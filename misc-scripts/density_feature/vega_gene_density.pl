#!/usr/local/bin/perl

=head1 NAME

vega_gene_density.pl -
script to calculate gene densities and stats in Vega

=head1 SYNOPSIS

    ./vega_gene_density.pl
        --species=Homo_sapiens
        [--dry-run|-n]
        [--help|-h]

=head1 DESCRIPTION

This script calculates Vega gene densities and total numbers per chromosome
for use in mapview.

=head1 LICENCE

This code is distributed under an Apache style licence:
Please see http://www.ensembl.org/code_licence.html for details

=head1 AUTHOR

Patrick Meidl <pm2@sanger.ac.uk>

based on code by
    Graham McVicer <mcvicker@ebi.ac.uk>
    Steve Trevanion <st3@sanger.ac.uk>

=head1 CONTACT

Post questions to the EnsEMBL development list ensembl-dev@ebi.ac.uk

=cut

use strict;
BEGIN {
    $ENV{'ENSEMBL_SERVERROOT'} = "../../..";
    unshift(@INC,"$ENV{'ENSEMBL_SERVERROOT'}/conf");
    unshift(@INC,"$ENV{'ENSEMBL_SERVERROOT'}/ensembl-compara/modules");
    unshift(@INC,"$ENV{'ENSEMBL_SERVERROOT'}/ensembl-draw/modules");
    unshift(@INC,"$ENV{'ENSEMBL_SERVERROOT'}/ensembl-external/modules");
    unshift(@INC,"$ENV{'ENSEMBL_SERVERROOT'}/ensembl-otter/modules");
    unshift(@INC,"$ENV{'ENSEMBL_SERVERROOT'}/modules");
    unshift(@INC,"$ENV{'ENSEMBL_SERVERROOT'}/ensembl/modules");
    unshift(@INC,"$ENV{'ENSEMBL_SERVERROOT'}/bioperl-live");
}

use SiteDefs;
use EnsWeb;
use EnsEMBL::DB::Core;
use Getopt::Long;
use Bio::EnsEMBL::DensityType;
use Bio::EnsEMBL::DensityFeature;
use POSIX;

use Data::Dumper;

my ($species, $dry, $help);
&GetOptions(
    "species=s" => \$species,
    "dry-run"   => \$dry,
    "n"         => \$dry,
    "help"      => \$help,
    "h"         => \$help,
);

if($help || !$species){
    print qq(Usage:
    ./vega_gene_density.pl
        --species=Homo_sapiens
        [--dry-run|-n]
        [--help|-h]\n\n);
    exit;
}

$ENV{'ENSEMBL_SPECIES'} = $species;

## set db user/pass to allow write access
my $db_ref = $EnsWeb::species_defs->databases;
$db_ref->{'ENSEMBL_DB'}{'USER'} = $EnsWeb::species_defs->ENSEMBL_WRITE_USER;
$db_ref->{'ENSEMBL_DB'}{'PASS'} = $EnsWeb::species_defs->ENSEMBL_WRITE_PASS;

## connect to databases
my $databases = &EnsEMBL::DB::Core::get_databases(qw(core));
my $db = $databases->{'core'};

die "Problem connecting to databases: $databases->{'error'}\n"
    if  $databases->{'error'} ;
warn "Database error: $databases->{'non_fatal_error'}\n"
    if $databases->{'non_fatal_error'};

## get adaptors
my $dfa = $db->get_DensityFeatureAdaptor;
my $dta = $db->get_DensityTypeAdaptor;
my $aa  = $db->get_AnalysisAdaptor;
my $attrib_adaptor = $db->get_AttributeAdaptor;
my $slice_adaptor = $db->get_SliceAdaptor;
my $top_slices = $slice_adaptor->fetch_all('toplevel');

## determine blocksize
my ($block_count, $genome_size, $block_size);
my $sth = $db->prepare( "select count(*) from gene" );
$sth->execute;
my ($gene_count) = $sth->fetchrow_array;
if( ! $gene_count ) {
    print STDERR "No gene density for " . $db->dbname . ".\n";
    exit;
} else {
    $block_count = $gene_count >> 1;
}
my @chr;
for my $slice ( @$top_slices ) {
    $genome_size += $slice->length;
    push @chr, $slice->seq_region_name;
}
$block_size = int( $genome_size / $block_count );
print STDERR "\nAvailable chromosomes: @chr\n";

## gene types
my %gene_types = (
    "Known" => "knownGeneDensity",
    "Novel_CDS" => "novelCDSDensity",
    "Novel_Transcript" => "novelTransDensity",
    "Putative" => "putativeTransDensity",
    "Predicted_Gene" => "predictedTransDensity",
    "Ig_Pseudogene_Segment" => "IgPseudoSegDensity",
    "Ig_Segment" => "IgSegDensity",
    "Processed_pseudogene" => "pseudoGeneDensity",
    "Pseudogene" => "pseudoGeneDensity",
    "Unprocessed_pseudogene" => "pseudoGeneDensity",
    "Known_in_progress" => "knownGeneDensity",
    "Novel_CDS_in_progress" => "novelCDSDensity", 
);
my %density_types;
map { $density_types{$_} = 1 } values %gene_types;

print STDERR "\nAvailable gene types: ";
print STDERR join(" ", sort keys %gene_types);
print STDERR "\n";

## create analysis and density type objects
my %dtcache;
foreach my $type (keys %gene_types) {
    my $analysis = new Bio::EnsEMBL::Analysis (
            -program     => "vega_gene_density_calc.pl",
            -database    => "ensembl",
            -gff_source  => "vega_gene_density_calc.pl",
            -gff_feature => "density",
            -logic_name  => $gene_types{$type});
    $aa->store($analysis) unless $dry;
    $analysis = $aa->fetch_by_logic_name($gene_types{$type});

    my $dt = Bio::EnsEMBL::DensityType->new(-analysis   => $analysis,
            -block_size => $block_size,
            -value_type => 'sum');
    $dta->store($dt) unless $dry;
    $dtcache{$gene_types{$type}} = $dt;
}

## loop over chromosomes
my ( $current_start, $current_end );
foreach my $slice (@$top_slices){
    $current_start = 1;
    my $chr = $slice->seq_region_name;
    my (%total, $i, %gene_names);
    my $bins = POSIX::ceil($slice->end / $block_size);
    
    print STDERR "\nGene densities for chr $chr with block size $block_size\n";
    print STDERR "Start at " . `date`;

    ## loop over blocks
    my @density_features;
    while($current_start <= $slice->end) {
        $i++;
        $current_end = $current_start+$block_size-1;
        if( $current_end > $slice->end ) {
            $current_end = $slice->end;
        }
        my $sub_slice = $slice->sub_Slice( $current_start, $current_end );
        my %num = ();

        ## count genes by type
        my $genes;
        eval { $genes = $sub_slice->get_all_Genes; };
        if ($@) {
            warn $@;
            $current_start = $current_end + 1;
            next;
        }
        foreach my $gene (@{$genes}) {
            ## only count genes that don't overlap the subslice start
            ## (since these were already counted in the last bin)
            if ($gene->start >= 1) {
	        $total{$gene->type}++;
            }
	    $num{$gene_types{$gene->type}}++;

        }

        ## create DensityFeature objects for each type
        foreach my $type (keys %density_types) {
  
            push @density_features, Bio::EnsEMBL::DensityFeature->new
                (-seq_region    => $slice,
                 -start         => $current_start,
                 -end           => $current_end,
                 -density_type  => $dtcache{$type},
                 -density_value => $num{$type} ||0
            );
        }
        $current_start = $current_end + 1;

        ## logging
        print STDERR "Chr: $chr | Bin: $i/$bins | Counts: ";
        print STDERR join(",", map { $num{$gene_types{$_}} || 0 }
                            sort keys %gene_types);
        print STDERR " | ";
        print STDERR "Mem: " . `ps $$ -o vsz |tail -1`;
    }


    ## store DensityFeatures for the chromosome
    $dfa->store(@density_features) unless $dry;

    ## stats
    my @attribs;
    push @attribs, Bio::EnsEMBL::Attribute->new
    (-NAME => 'Known genes',
     -CODE => 'KnownGeneCount',
     -VALUE => $total{'Known'} || 0,
     -DESCRIPTION => 'Total Number of Known genes');

    push @attribs, Bio::EnsEMBL::Attribute->new
    (-NAME => 'Novel CDS',
     -CODE => 'NovelCDSCount',
     -VALUE => $total{'Novel_CDS'} || 0,
     -DESCRIPTION => 'Total Number of Novel CDSs');

    push @attribs, Bio::EnsEMBL::Attribute->new
    (-NAME => 'Novel transcripts',
     -CODE => 'NovelTransCount',
     -VALUE => $total{'Novel_Transcript'} || 0,
     -DESCRIPTION => 'Total Number of Novel transcripts');

    push @attribs, Bio::EnsEMBL::Attribute->new
    (-NAME => 'Putative transcripts',
     -CODE => 'PutTransCount',
     -VALUE => $total{'Putative'} || 0,
     -DESCRIPTION => 'Total Number of Putative transcripts');

    push @attribs, Bio::EnsEMBL::Attribute->new
    (-NAME => 'Predicted transcripts',
     -CODE => 'PredTransCount',
     -VALUE => $total{'Predicted_Gene'} || 0,
     -DESCRIPTION => 'Total Number of Predicted transcripts');

    push @attribs, Bio::EnsEMBL::Attribute->new
    (-NAME => 'Ig Segments',
     -CODE => 'IgSegCount',
     -VALUE => $total{'Ig_Segment'} || 0,
     -DESCRIPTION => 'Total Number of Ig Segments');

    push @attribs, Bio::EnsEMBL::Attribute->new
    (-NAME => 'Ig Pseudogene Segments',
     -CODE => 'IgPsSegCount',
     -VALUE => $total{'Ig_Pseudogene_Segment'} || 0,
     -DESCRIPTION => 'Total Number of Ig Pseudogene Segments');

    push @attribs, Bio::EnsEMBL::Attribute->new
    (-NAME => 'Total Pseudogenes',
     -CODE => 'TotPsCount',
     -VALUE => $total{'Pseudogenes'}
               + $total{'Processed_pseudogene'}
               + $total{'Unprocessed_pseudogene' || 0},
     -DESCRIPTION => 'Total Number of Pseudogenes');

    push @attribs, Bio::EnsEMBL::Attribute->new
    (-NAME => 'Unclassified pseudogenes',
     -CODE => 'UnclassPsCount',
     -VALUE => $total{'Pseudogene'} || 0,
     -DESCRIPTION => 'Number of Unclassified pseudogenes');

    push @attribs, Bio::EnsEMBL::Attribute->new
    (-NAME => 'Processed pseudogenes',
     -CODE => 'ProcPsCount',
     -VALUE => $total{'Processed_pseudogene'} || 0,
     -DESCRIPTION => 'Number of Processed pseudogenes');

    push @attribs, Bio::EnsEMBL::Attribute->new
    (-NAME => 'Unprocessed pseudogenes',
     -CODE => 'UnprocPsCount',
     -VALUE => $total{'Unprocessed_pseudogene'} || 0,
     -DESCRIPTION => 'Number of Unprocessed pseudogenes');

#    $attrib_adaptor->store_on_Slice($slice, \@attribs) unless $dry;

    print STDERR "Total for chr $chr:\n";
    print STDERR map { "\t$_ => $total{$_}\n" } sort keys %total;

}

print STDERR "\nAll done at " . `date` . "\n";

