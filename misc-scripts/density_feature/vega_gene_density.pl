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

my ($species, $chr, $dry, $help);
&GetOptions(
    "species=s" => \$species,
    "chr=s"     => \$chr,
    "dry-run"   => \$dry,
    "n"         => \$dry,
    "help"      => \$help,
    "h"         => \$help,
);

if($help || !$species){
    print qq(Usage:
    ./vega_gene_density.pl
        --species=Homo_sapiens
        [--chr=1,2]
        [--dry-run|-n]
        [--help|-h]\n\n);
    exit;
}

$ENV{'ENSEMBL_SPECIES'} = $species;

#get the adaptors needed
my $slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species,"vega","Slice") or die "can't load slice adaptor - is the species name correct?";
my $dfa =  Bio::EnsEMBL::Registry->get_adaptor($species,"vega","DensityFeature") or die;
my $dta =  Bio::EnsEMBL::Registry->get_adaptor($species,"vega","DensityType") or die;
my $aa  =  Bio::EnsEMBL::Registry->get_adaptor($species,"vega","Analysis") or die;
my $attrib_adaptor =  Bio::EnsEMBL::Registry->get_adaptor($species,"vega","Attribute");

## set db user/pass to allow write access
$EnsWeb::species_defs->set_write_access('ENSEMBL_DB',$species);

## which chromosomes do we run?
my @top_slices;
if ($chr) {
    ## run chromosomes specified on commandline
    foreach (split(",", $chr)) {
        push @top_slices, $slice_adaptor->fetch_by_region("toplevel", $_);
    }
} else {
    ## run all chromosomes for this species
    @top_slices = @{$slice_adaptor->fetch_all("toplevel")};
}

## determine blocksize, assuming you want 150 blocks for the smallest
## chromosome over 5Mb in size. Use an extra, smaller bin size for chromosomes smaller than 5Mb
my $big_chr = [];
my $small_chr = [];
my (@big_chr_names, $big_block_size, $min_big_chr);
my (@small_chr_names, $small_block_size, $min_small_chr);
for my $slice ( @top_slices ) {
    if ($slice->length < 5000000) {
	if (! $min_small_chr or ($min_small_chr > $slice->length)) {
	    $min_small_chr = $slice->length;
	}
	push @small_chr_names, $slice->seq_region_name;
	push @{$small_chr->[0]}, $slice;
    }
    if (! $min_big_chr or ($min_big_chr > $slice->length) && $slice->length > 5000000) {
	$min_big_chr = $slice->length;
    }
    push @big_chr_names, $slice->seq_region_name;
    push @{$big_chr->[0]}, $slice;
}

$big_block_size = int( $min_big_chr / 150 );
push @{$big_chr}, $big_block_size;
$small_block_size = int( $min_small_chr / 150 );
push @{$small_chr}, $small_block_size;

print STDERR "\nAvailable chromosomes using block size of $big_block_size: @big_chr_names\n";
print STDERR "\nAvailable chromosomes using block size of $small_block_size: @small_chr_names\n";

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
foreach my $block_size ($big_block_size,$small_block_size) {
    eval {
	foreach my $type (keys %gene_types) {
	    my $analysis = new Bio::EnsEMBL::Analysis (
						       -program     => "vega_gene_density.pl",
						       -database    => "ensembl",
						       -gff_source  => "vega_gene_density.pl",
						       -gff_feature => "density",
						       -logic_name  => $gene_types{$type});
	    $aa->store($analysis) unless $dry;
	    $analysis = $aa->fetch_by_logic_name($gene_types{$type});
	    my $dt = Bio::EnsEMBL::DensityType->new(-analysis   => $analysis,
						    -block_size => $block_size,
						    -value_type => 'sum');
	    $dta->store($dt) unless $dry;
	    $dtcache{$block_size}{$gene_types{$type}} = $dt;
	}    
    }
}


## loop over chromosomes, doing big ones then small ones
my ( $current_start, $current_end );
foreach my $object ($big_chr, $small_chr) {
    eval {
	my $block_size = $object->[1];
	foreach my $slice (@{$object->[0]}){
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
			 -density_type  => $dtcache{$block_size}{$type},
			 -density_value => $num{$type} ||0
			);
		}
		$current_start = $current_end + 1;
		
		## logging
		print STDERR "Chr: $chr | Bin: $i/$bins | Counts: ";
		print STDERR join(",", map { $num{$gene_types{$_}} || 0 }
				  sort keys %gene_types);
		print STDERR " | ";
		print STDERR "Mem: " . `ps -p $$ -o vsz |tail -1`;
	    }
	    
	    
	    ## store DensityFeatures for the chromosome
	    $dfa->store(@density_features) unless $dry;
	    
	    ## stats
	    my @attribs;
	    push @attribs, Bio::EnsEMBL::Attribute->new
		(-NAME => '12:Known genes',
		 -CODE => 'KnownGeneCount',
		 -VALUE => $total{'Known'} || 0,
		 -DESCRIPTION => 'Total Number of Known genes');
	    
	    push @attribs, Bio::EnsEMBL::Attribute->new
		(-NAME => '14:Novel CDS',
		 -CODE => 'NovelCDSCount',
		 -VALUE => $total{'Novel_CDS'} || 0,
		 -DESCRIPTION => 'Total Number of Novel CDSs');
	    
	    push @attribs, Bio::EnsEMBL::Attribute->new
		(-NAME => '16:Novel transcripts',
		 -CODE => 'NovelTransCount',
		 -VALUE => $total{'Novel_Transcript'} || 0,
		 -DESCRIPTION => 'Total Number of Novel transcripts');
	    
	    push @attribs, Bio::EnsEMBL::Attribute->new
		(-NAME => '24:Putative',
		 -CODE => 'PutTransCount',
		 -VALUE => $total{'Putative'} || 0,
		 -DESCRIPTION => 'Total Number of Putative transcripts');
	    
	    push @attribs, Bio::EnsEMBL::Attribute->new
		(-NAME => '30:Predicted transcripts',
		 -CODE => 'PredTransCount',
		 -VALUE => $total{'Predicted_Gene'} || 0,
		 -DESCRIPTION => 'Total Number of Predicted transcripts');
	    
	    push @attribs, Bio::EnsEMBL::Attribute->new
		(-NAME => '26:Ig segments',
		 -CODE => 'IgSegCount',
		 -VALUE => $total{'Ig_Segment'} || 0,
		 -DESCRIPTION => 'Total Number of Ig Segments');
	    
	    push @attribs, Bio::EnsEMBL::Attribute->new
		(-NAME => '28:Ig pseudogene segments',
		 -CODE => 'IgPsSegCount',
		 -VALUE => $total{'Ig_Pseudogene_Segment'} || 0,
		 -DESCRIPTION => 'Total Number of Ig Pseudogene Segments');
	    
	    push @attribs, Bio::EnsEMBL::Attribute->new
		(-NAME => '18:Total pseudogenes',
		 -CODE => 'TotPsCount',
		 -VALUE => $total{'Pseudogene'}
		 + $total{'Processed_pseudogene'}
		 + $total{'Unprocessed_pseudogene' || 0},
		 -DESCRIPTION => 'Total Number of Pseudogenes');
	    
	    push @attribs, Bio::EnsEMBL::Attribute->new
		(-NAME => '20:Processed pseudogenes',
		 -CODE => 'ProcPsCount',
		 -VALUE => $total{'Processed_pseudogene'} || 0,
		 -DESCRIPTION => 'Number of Processed pseudogenes');
	    
	    push @attribs, Bio::EnsEMBL::Attribute->new
		(-NAME => '22:Unprocessed pseudogenes',
		 -CODE => 'UnprocPsCount',
		 -VALUE => $total{'Unprocessed_pseudogene'} || 0,
		 -DESCRIPTION => 'Number of Unprocessed pseudogenes');
	    
	    push @attribs, Bio::EnsEMBL::Attribute->new
		(-NAME => '13:Known genes (in progress)',
		 -CODE => 'KnwnprogCount',
		 -VALUE => $total{'Known_in_progress'} || 0,
		 -DESCRIPTION => 'Number of Known Genes in progress');
	    
	    push @attribs, Bio::EnsEMBL::Attribute->new
		(-NAME => '15:Novel CDS (in progress)',
		 -CODE => 'NovCDSprogCount',
		 -VALUE => $total{'Novel_CDS_in_progress'} || 0,
		 -DESCRIPTION => 'Number of novel CDS in progress');
	    
	    #only store unclassified pseudogenes if there are no processed and unprocessed pseudos, ie if
	    #total pseudos eq pseudos
	    unless ($total{'Unprocessed_pseudogene'} == 0 && $total{'Processed_pseudogene'} == 0) {  
		push @attribs, Bio::EnsEMBL::Attribute->new
		    (-NAME => '23:Unclassified pseudogenes',
		     -CODE => 'UnclassPsCount',
		     -VALUE => $total{'Pseudogene'} || 0,
		     -DESCRIPTION => 'Number of Unclassified pseudogenes');
	    }
	    
	    $attrib_adaptor->store_on_Slice($slice, \@attribs) unless $dry;
	    
	    print STDERR "Total for chr $chr:\n";
	    print STDERR map { "\t$_ => $total{$_}\n" } sort keys %total;
	}
    }
}
print STDERR "\nAll done at " . `date` . "\n";

