#!/usr/local/ensembl/bin/perl -w
#
# Calculate the repeat coverage for given database.
# condition: 1k blocks to show contigview displays
#  4000 blocks for a whole genome
#
# checks wether database contains repeats before doing anything

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::DBSQL::DensityFeatureAdaptor;
use Bio::EnsEMBL::DBSQL::DensityTypeAdaptor;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::DensityType;
use Bio::EnsEMBL::DensityFeature;
use Bio::EnsEMBL::Mapper::RangeRegistry;
use Bio::EnsEMBL::Utils::Exception qw(warning throw);

use POSIX;

use strict;
use Getopt::Long;

my ( $host, $user, $pass, $port, $dbname  );

GetOptions( "host=s", \$host,
	    "user=s", \$user,
	    "pass=s", \$pass,
	    "port=i", \$port,
	    "dbname=s", \$dbname
	  );


my $chunksize = 1_000_000;
my $small_blocksize = 1_000;
my $bin_count = 150;
my $max_top_slice = 100;

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $host,
					    -user => $user,
					    -port => $port,
					    -pass => $pass,
					    -dbname => $dbname);



my $sth = $db->dbc()->prepare( "select count(*) from repeat_feature" );
$sth->execute();

my ( $repeat_count )  = $sth->fetchrow_array();

if( ! $repeat_count ) {
  print STDERR "No repeat density for $dbname.\n";
  exit();
}

#
# Get the adaptors needed;
#

#
# Clean up old features first. Also remove analysis and density type entry as these are recreated.
#

print "Deleting old PercentageRepeat features\n";
$sth = $db->dbc->prepare("DELETE df, dt, a, ad FROM analysis_description ad, density_feature df, density_type dt, analysis a WHERE ad.analysis_id = a.analysis_id AND a.analysis_id=dt.analysis_id AND dt.density_type_id=df.density_type_id AND a.logic_name='PercentageRepeat'");
$sth->execute();

# $sth = $db->dbc()->prepare(
#   qq(
#   DELETE ad
#   FROM analysis_description ad
#   WHERE ad.display_label = 'PercentageRepeat') );
# $sth->execute();

my $slice_adaptor = $db->get_SliceAdaptor();
my $dfa = $db->get_DensityFeatureAdaptor();
my $dta = $db->get_DensityTypeAdaptor();
my $aa  = $db->get_AnalysisAdaptor();



#
# Create new analysis object for density calculation.
#

my $analysis =
  new Bio::EnsEMBL::Analysis(
       -program     => "repeat_coverage_calc.pl",
       -database    => "ensembl",
       -gff_source  => "repeat_coverage_calc.pl",
       -gff_feature => "density",
       -logic_name  => "PercentageRepeat",
       -description =>
         'Percentage of repetetive elements for top level seq_regions.',
       -display_label => 'Repeats (percent)',
       -displayable   => 1 );

$aa->store($analysis);
$aa->update($analysis);

my $slices = $slice_adaptor->fetch_all( "toplevel" );
my @sorted_slices = sort { $b->seq_region_length() <=> $a->seq_region_length() } @$slices;

my $small_density_type = Bio::EnsEMBL::DensityType->new
  (-analysis   => $analysis,
   -block_size => $small_blocksize,
   -value_type => 'ratio');

my $variable_density_type = Bio::EnsEMBL::DensityType->new
  (-analysis   => $analysis,
   -region_features => $bin_count,
   -value_type => 'ratio');

$dta->store($small_density_type);
$dta->store($variable_density_type);


my $slice_count = 0;



foreach my $slice ( @sorted_slices ) {
  
  #
  # do it for small and large blocks
  #
  print STDERR ("Working on seq_region ".$slice->seq_region_name()." length ".$slice->seq_region_length());

  my $rr = Bio::EnsEMBL::Mapper::RangeRegistry->new();
  my $chunk_end = 0;
  my $variable_end = 0;
  my $small_end = 0;
  my ( $small_start );
  my $repeat_size;
  my $variable_start = 0;
  my $variable_blocksize = POSIX::ceil( $slice->seq_region_length() / 
					$variable_density_type->region_features());
  $slice_count++;

  while( $chunk_end < $slice->length() ) {
    my $chunk_start = $chunk_end+1;
    $chunk_end += $chunksize;
    $chunk_end = $slice->length() if $chunk_end > $slice->length();

    register( $rr, $slice, $chunk_start, $chunk_end );

    my @dfs = ();

    if( $slice_count < $max_top_slice ) {
      while ( $variable_end+$variable_blocksize <= $chunk_end ) {
	# here we can do the variable sized repeat densities
	$variable_start = $variable_end+1;
	$variable_end += $variable_blocksize;

	$repeat_size = $rr->overlap_size( "1", $variable_start, $variable_end );
	my $percentage_repeat = $repeat_size / $variable_blocksize * 100;

	push( @dfs, Bio::EnsEMBL::DensityFeature->new
	      (-seq_region    => $slice,
	       -start         => $variable_start,
	       -end           => $variable_end,
	       -density_type  => $variable_density_type,
	       -density_value => $percentage_repeat));

      }
    }

    while ( $small_end + $small_blocksize <= $chunk_end ) {
      # here we can do the small sized density features
      $small_start = $small_end+1;
      $small_end += $small_blocksize;

      $repeat_size = $rr->overlap_size( "1", $small_start, $small_end );
      my $percentage_repeat = $repeat_size / $small_blocksize * 100;

      push( @dfs, Bio::EnsEMBL::DensityFeature->new
        (-seq_region    => $slice,
         -start         => $small_start,
         -end           => $small_end,
         -density_type  => $small_density_type,
         -density_value => $percentage_repeat));
    }

    if (@dfs) {
      $dfa->store(@dfs);
    } else {
      warning("No repeat density calculated for ".$slice->name." (chunk start $chunk_start, chunk end $chunk_end).");
    }

    my $used_lower_limit = $small_start<$variable_start?$small_start:$variable_start;

    # here some rr cleanup
    $rr->check_and_register( "1", 0, $used_lower_limit );
  }

  # missing the last bits
  if( $small_end < $slice->length() ) {
    $small_start = $small_end+1;
    $small_end = $slice->length();

    $repeat_size = $rr->overlap_size( "1", $small_start, $small_end );
    my $percentage_repeat = $repeat_size / ($small_end - $small_start + 1 ) * 100;

    $dfa->store( Bio::EnsEMBL::DensityFeature->new
		 (-seq_region    => $slice,
		  -start         => $small_start,
		  -end           => $small_end,
		  -density_type  => $small_density_type,
		  -density_value => $percentage_repeat));
  }

  if( $variable_end < $slice->length() && $slice_count < $max_top_slice ) {
    $variable_start = $variable_end+1;
    $variable_end = $slice->length();

    $repeat_size = $rr->overlap_size( "1", $variable_start, $variable_end );
    my $percentage_repeat = $repeat_size / ($variable_end - $variable_start+1) * 100;

    $dfa->store( Bio::EnsEMBL::DensityFeature->new
		 (-seq_region    => $slice,
		  -start         => $variable_start,
		  -end           => $variable_end,
		  -density_type  => $variable_density_type,
		  -density_value => $percentage_repeat));
  }

  print STDERR " DONE.\n";
}


sub register {
  my ($rr, $slice, $start, $end ) = @_;

  my $subSlice = $slice->sub_Slice( $start, $end );
  my $repeats = $subSlice->get_all_RepeatFeatures();
  for my $repeat ( @$repeats ) {
    $rr->check_and_register( "1", $repeat->seq_region_start(), $repeat->seq_region_end() );
  }
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
#  my $avg = undef;
#  $avg = $sum/$length if($length < 0);
#  print("Sum=$sum, Length=$length, Avg/Base=$sum");
}




  


