#!/usr/local/bin/perl

=head1 NAME

glovar_snp_density.pl -
Script to calculate glovar SNP density and stats (and optioanlly prepare AV
index dumps) for Vega.

=head1 SYNOPSIS

    ./glovar_snp_density.pl
        --species=Homo_sapiens
        [--chr=6,13,14]
        [--dry-run|-n]
        [--avdump|-a]
        [--help|-h]

=head1 DESCRIPTION

This script calculates Glovar SNP densities and total numbers per chromosome
for use in mapview. Can be run for individual chromosomes if desired (default:
all chromosomes). It optionally also dumps SNPs into a file for generating the
AV search index.

=head1 LICENCE

This code is distributed under an Apache style licence:
Please see http://www.ensembl.org/code_licence.html for details

=head1 AUTHOR

Patrick Meidl <pm2@sanger.ac.uk>

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
    unshift(@INC,"$ENV{'ENSEMBL_SERVERROOT'}/ensembl-variation/modules");
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
use Bio::EnsEMBL::Registry;
my $reg = "Bio::EnsEMBL::Registry";

my ($species, $chr, $dry, $avdump, $help);
&GetOptions(
    "species=s" => \$species,
    "chr=s"     => \$chr,
    "dry-run"   => \$dry,
    "n"         => \$dry,
    "avdump"    => \$avdump,
    "a"         => \$avdump,
    "help"      => \$help,
    "h"         => \$help,
);

if($help || !$species){
    print qq(Usage:
    ./glovar_snp_density.pl
        --species=Homo_sapiens
        [--chr=6,13,14]
        [--dry-run|-n]
        [--avdump|-a]
        [--help|-h]\n\n);
    exit;
}

$ENV{'ENSEMBL_SPECIES'} = $species;

## set db user/pass to allow write access
$EnsWeb::species_defs->set_write_access('ENSEMBL_DB',$species);

# connect to databases
my $databases = &EnsEMBL::DB::Core::get_databases(qw(core glovar));

die "Problem connecting to databases: $databases->{'error'}\n"
    if  $databases->{'error'} ;
warn "Database error: $databases->{'non_fatal_error'}\n"
    if $databases->{'non_fatal_error'};

# get the adaptors needed
my $dfa = $databases->{'core'}->get_DensityFeatureAdaptor;
my $dta = $databases->{'core'}->get_DensityTypeAdaptor;
my $aa  = $databases->{'core'}->get_AnalysisAdaptor;
my $attrib_adaptor = $databases->{'core'}->get_AttributeAdaptor;
my $slice_adaptor = $databases->{'core'}->get_SliceAdaptor;

# which chromosomes do we run?
my @top_slices;
if ($chr) {
    # run chromosomes specified on commandline
    foreach (split(",", $chr)) {
        push @top_slices, $slice_adaptor->fetch_by_region("toplevel", $_);
    }
} else {
    # run all chromosomes for this species
    @top_slices = @{$slice_adaptor->fetch_all("toplevel")};
}

# calculate block size (assuming 4000 blocks per genome)
my ( $block_size, $genome_size );
for my $slice ( @{$slice_adaptor->fetch_all("toplevel")} ) {
    $genome_size += $slice->length;
}
$block_size = int( $genome_size / 4000 );

# analysis
my $analysis = new Bio::EnsEMBL::Analysis (
        -program     => "glovar_snp_density.pl",
        -database    => "vega",
        -gff_source  => "glovar_snp_density.pl",
        -gff_feature => "density",
        -logic_name  => "snpDensity");
$aa->store( $analysis ) unless $dry;

# density type
my $dt = Bio::EnsEMBL::DensityType->new(
        -analysis   => $analysis,
        -block_size => $block_size,
        -value_type => 'sum');
$dta->store($dt) unless $dry;

# loop over chromosomes
my @chr;
foreach my $sl (@top_slices) {
    push @chr, $sl->seq_region_name;
}
print STDERR "\nAvailable chromosomes: @chr\n";

# settings for AV index dump
use constant SNP_LINE => join("\t",
    'Glovar SNP', '%s', '/%s/snpview?snp=%s&source=glovar', '%s',
    "Single nucleotide polymorphism (SNP) %s [Alleles: %s]. Alternative IDs: %s.\n"
);
if ($avdump) {
    print STDERR "Preparing directories for AV index dump...\n";
    my $dumpdir = "$ENV{'ENSEMBL_SERVERROOT'}/utils/indexing/input";
    unless (-e $dumpdir) {
        mkdir $dumpdir, 0777 or die "Could not creat directory $dumpdir: $!\n";
    }
    unless (-e "$dumpdir/$species") {
        mkdir "$dumpdir/$species", 0777 or die
            "Could not creat directory $dumpdir/$species: $!\n";
    }
    open (AV, ">>$dumpdir/$species/SNP.txt") or die
        "Could not open $dumpdir/$species/SNP.txt for writing: $!\n";
    print STDERR "Done.\n";
}

my ($current_start, $current_end);
foreach my $slice (@top_slices) {
    $current_start = 1;
    my $chr = $slice->seq_region_name;
    my ($total, $i);
    my $bins = POSIX::ceil($slice->end / $block_size);
    
    print STDERR "\nSNP densities for chr $chr with block size $block_size\n";
    print STDERR "Start at " . `date`;

    # loop over blocks
    while($current_start <= $slice->end) {
        $i++;
        $current_end = $current_start+$block_size-1;
        if ($current_end > $slice->end) {
            $current_end = $slice->end;
        }
        my $sub_slice = $slice->sub_Slice( $current_start, $current_end );
        my $count = 0;

        my $varfeats;
        eval { $varfeats = $sub_slice->get_all_ExternalFeatures('GlovarSNP'); };
        if ($@) {
            warn $@;
            $current_start = $current_end + 1;
            next;
        }
        # only count varfeats that don't overlap slice start
        # also, avoid duplicate counting
        my %varfeats = map { "$_->variation_name => 1" if ($_->start >= 1) } @{$varfeats};
        $count = scalar(keys %varfeats);

        # AV index dump
        if ($avdump) {
            foreach my $varfeat (@{$varfeats}) {
                next if ($varfeat->start < 1);
                my $snpid = $varfeat->variation_name;

                # dblinks
                my @sources = @{ $varfeat->variation->get_all_synonym_sources };
                my (@IDs, @desc);
                foreach my $source (@sources) {
                    my @extIDs = @{ $varfeat->variation->get_all_synonyms($source) };
                    push @IDs, @extIDs;
                    push @desc, "$source: @extIDs";
                }
                
                print AV sprintf SNP_LINE,
                    $snpid,
                    $species,
                    $snpid,
                    join(" ", @IDs),
                    $snpid,
                    $varfeat->allele_string,
                    join(", ", @desc)
                ; 
            }
        }

        # density
        my $df = Bio::EnsEMBL::DensityFeature->new
            (-seq_region    => $slice,
             -start         => $current_start,
             -end           => $current_end,
             -density_type  => $dt,
             -density_value => $count);
        $current_start = $current_end + 1;
        $dfa->store($df) unless $dry;
        $total += $count;

        # logging
        print STDERR "Chr: $chr | Bin: $i/$bins | Count: $count | ";
        print STDERR "Mem: " . `ps -p $$ -o vsz |tail -1`;
    }

    # stats
    print STDERR "Total for chr $chr: $total\n";
    my $stat = Bio::EnsEMBL::Attribute->new
	(-NAME => 'SNPs',
	 -CODE => 'SNPCount',
	 -VALUE => $total,
	 -DESCRIPTION => 'Total Number of SNPs');
    $attrib_adaptor->store_on_Slice($slice, [$stat]) unless $dry;
}
close AV if $avdump;

print STDERR "\nAll done at " . `date` . "\n";

