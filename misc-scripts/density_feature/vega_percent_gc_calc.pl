#!/usr/bin/perl -w

#
# Calculate the GC content for top level seq_regions
#   small regions 500bp to be able to display on contigview
#   big regions genomesize / 4000 for 4000 features on the genome

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

#
# Get the adaptors needed;
#
my $slice_adaptor = $db->get_SliceAdaptor();
my $dfa = $db->get_DensityFeatureAdaptor();
my $dta = $db->get_DensityTypeAdaptor();
my $aa  = $db->get_AnalysisAdaptor();

my $top_slices = $slice_adaptor->fetch_all( "toplevel" );

my $big_chr = [];
my $small_chr = [];

my (@big_chr_names, $big_block_size, $min_big_chr);
my (@small_chr_names, $small_block_size, $min_small_chr);
for my $slice ( @$top_slices ) {
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
$small_block_size = int( $min_small_chr / 150 );
push @{$big_chr}, $big_block_size;
push @{$small_chr}, $small_block_size;


print STDERR "\nAvailable chromosomes using block size of $big_block_size: @big_chr_names\n";
print STDERR "\nAvailable chromosomes using block size of $small_block_size: @small_chr_names\n";

#
# Create new analysis object for density calculation.
#
my $analysis = new Bio::EnsEMBL::Analysis (-program     => "percent_gc_calc.pl",
					   -database    => "ensembl",
					   -gff_source  => "percent_gc_calc.pl",
					   -gff_feature => "density",
					   -logic_name  => "PercentGC");
$aa->store($analysis) unless $dry;

#
# Create new density types
#
if ($small_block_size) {
    my $small_density_type = Bio::EnsEMBL::DensityType->new
	(-analysis   => $analysis,
	 -block_size => $small_block_size,
	 -value_type => 'ratio');
    $dta->store($small_density_type) unless $dry;
    push @{$small_chr}, $small_density_type;
}

if ($big_block_size) {
    my $big_density_type = Bio::EnsEMBL::DensityType->new
	(-analysis   => $analysis,
	 -block_size => $big_block_size,
	 -value_type => 'ratio');
    $dta->store($big_density_type) unless $dry;
    push @{$big_chr}, $big_density_type;
}

my ( $current_start, $current_end );

$Data::Dumper::Maxdepth = 2;

my $types = [];
foreach my $object ( $big_chr, $small_chr) {
    eval {
	my $block_size = $object->[1];
	foreach my $slice (@{$object->[0]}){
	    $current_start = 1;
	    warn "Chromosome ",$slice->seq_region_name;
	    while($current_start <= $slice->end()) {
		$current_end = $current_start+$block_size-1;
		if( $current_end > $slice->end() ) {
		    $current_end = $slice->end();
		}
		my $sub_slice = $slice->sub_Slice( $current_start, $current_end );
		warn "start = $current_start, end = $current_end\n";
		my $gc = $sub_slice->get_base_count()->{'%gc'};
		my $df =  Bio::EnsEMBL::DensityFeature->new
		    (-seq_region    => $slice,
		     -start         => $current_start,
		     -end           => $current_end,
		     -density_type  => $object->[2],
		     -density_value => $gc);
		$dfa->store($df) unless $dry;
		$current_start = $current_end+1;
	    }
	}
    };
}

#
# helper to draw an ascii representation of the density features
#
sub print_features {
    my $features = shift;
  
    return if(!@$features);
    
    my $sum = 0;
    my $length = 0;
    #  my $type = $features->[0]->{'density_type'}->value_type();
    
    print("\n");
    my $max=0;
    foreach my $f (@$features) {
	if($f->density_value() > $max){
	    $max=$f->density_value();
	}
    }
    foreach my $f (@$features) {
	my $i=1;
	for(; $i< ($f->density_value()/$max)*40; $i++){
	    print "*";
	}
	for(my $j=$i;$j<40;$j++){
	    print " ";
	}
	print "  ".$f->density_value()."\t".$f->start()."\n";
    }
}
