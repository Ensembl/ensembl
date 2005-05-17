#!/usr/local/bin/perl -w
#
# Calculate the repeat coverage for given database.
# condition: 1k blocks to show contigview displays
#  4000 blocks for a whole genome
#
# checks wether database contains repeats before doing anything
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
use Bio::EnsEMBL::Mapper::RangeRegistry;
use POSIX;

use Data::Dumper;

my ($species, $dry, $help);
&GetOptions(
    "species=s" => \$species,
    "dry_run"   => \$dry,
    "n"         => \$dry,
    "help"      => \$help,
    "h"         => \$help,
);

if($help || !$species){
    print qq(Usage:
    ./vega_gene_density.pl
        --species=Homo_sapiens
        [--dry_run|-n]
        [--help|-h]\n\n);
    exit;
}

$ENV{'ENSEMBL_SPECIES'} = $species;

#get the adaptors needed
my $slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species,"vega","Slice") or die "can't load slice adaptor - is the species name correct?";
my $dfa =  Bio::EnsEMBL::Registry->get_adaptor($species,"vega","DensityFeature") or die;
my $dta =  Bio::EnsEMBL::Registry->get_adaptor($species,"vega","DensityType") or die;
my $aa  =  Bio::EnsEMBL::Registry->get_adaptor($species,"vega","Analysis") or die;

## set db user/pass to allow write access
$EnsWeb::species_defs->set_write_access('ENSEMBL_DB',$species);

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
my $analysis = new Bio::EnsEMBL::Analysis (-program     => "repeat_coverage_calc.pl",
					   -database    => "ensembl",
					   -gff_source  => "repeat_coverage_calc.pl",
					   -gff_feature => "density",
					   -logic_name  => "PercentageRepeat");
$aa->store($analysis) unless $dry;

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

foreach my $object ($big_chr, $small_chr) {
    eval {
	my $block_size = $object->[1];
	foreach my $slice ( @{$object->[0]} ) {
	    $current_start = 1;
	    warn "Chromosome ", $slice->seq_region_name;
	    while($current_start <= $slice->end()) {
		$current_end = $current_start+$block_size-1;
		if( $current_end > $slice->end() ) {
		    $current_end = $slice->end();
		}
		my $this_block_size = $current_end - $current_start + 1;
		my $sub_slice = $slice->sub_Slice( $current_start, $current_end );
		warn "start = $current_start, end = $current_end\n";
		my $rr = Bio::EnsEMBL::Mapper::RangeRegistry->new() or die;
		foreach my $repeat (@{$sub_slice->get_all_RepeatFeatures()}){
		    $rr->check_and_register("1",$repeat->start,$repeat->end);
		}
		my $count = 0;
		my $non_repeats = $rr->check_and_register("1",1,$this_block_size);
		if( defined $non_repeats ) {
		    foreach my $non_repeat ( @$non_repeats ) {
			$count += ($non_repeat->[1]-$non_repeat->[0])+1;
		    }
		}
		my $percentage_repeat = (($this_block_size-$count)/$this_block_size)*100;
		my $df = Bio::EnsEMBL::DensityFeature->new
		    (-seq_region    => $slice,
		     -start         => $current_start,
		     -end           => $current_end,
		     -density_type  => $object->[2],
		     -density_value => $percentage_repeat);
		$dfa->store($df) unless $dry;
		$current_start = $current_end + 1;
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
