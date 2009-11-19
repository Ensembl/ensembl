#!/usr/local/ensembl/bin/perl -w

#
# Calculate the GC content for top level seq_regions
#   small regions 500bp to be able to display on contigview
#   big regions genomesize / 4000 for 4000 features on the genome


use strict;
#use dbi;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::DBSQL::DensityFeatureAdaptor;
use Bio::EnsEMBL::DBSQL::DensityTypeAdaptor;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::DensityType;
use Bio::EnsEMBL::DensityFeature;
use Getopt::Long;

my ( $host, $user, $pass, $port, $dbname  );

GetOptions( "host=s", \$host,
	    "user=s", \$user,
	    "pass=s", \$pass,
	    "port=i", \$port,
	    "dbname=s", \$dbname
	  );



my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $host,
					    -user => $user,
					    -port => $port,
					    -pass => $pass,
					    -dbname => $dbname);


my $bin_count  = 150;
my $max_slices = 100;

#
# Check wether the script should run on given database
#
my $sth = $db->dbc->prepare( "select count(*) from dna" );
$sth->execute();
my ( $dna_count )  = $sth->fetchrow_array();
if( ! $dna_count ) {
  print STDERR "No dna, no gc content for $dbname.\n";
  exit();
}



print "Deleting old PercentGC features\n";
$sth = $db->dbc->prepare(
  qq(
DELETE df, dt, a, ad
FROM density_feature df, density_type dt, analysis a, analysis_description ad
WHERE a.analysis_id=dt.analysis_id
AND ad.analysis_id = a.analysis_id
AND dt.density_type_id=df.density_type_id
AND a.logic_name='PercentGC') );
$sth->execute();

# $sth = $db->dbc()->prepare(
#   qq(
#   DELETE ad
#   FROM analysis_description ad
#   WHERE ad.display_label = 'PercentGC') );
# $sth->execute();

#
# Get the adaptors needed;
#

my $slice_adaptor = $db->get_SliceAdaptor();
my $dfa = $db->get_DensityFeatureAdaptor();
my $dta = $db->get_DensityTypeAdaptor();
my $aa  = $db->get_AnalysisAdaptor();

# Sort slices by coordinate system rank, then by length.
my @sorted_slices =
  sort( {      $a->coord_system()->rank() <=> $b->coord_system()->rank()
            || $b->seq_region_length() <=> $a->seq_region_length()
  } @{ $slice_adaptor->fetch_all('toplevel') } );

#
# Create new analysis object for density calculation.
#

my $analysis =
  new Bio::EnsEMBL::Analysis(
             -program     => "percent_gc_calc.pl",
             -database    => "ensembl",
             -gff_source  => "percent_gc_calc.pl",
             -gff_feature => "density",
             -logic_name  => "PercentGC",
             -description => 'Percentage of G/C bases in the sequence.',
             -display_label => 'PercentGC',
             -displayable   => 1 );

$aa->store($analysis);
$aa->update($analysis);

#
# Create new density type.
#


my $density_type = Bio::EnsEMBL::DensityType->new
  (-analysis   => $analysis,
   -region_features => $bin_count,
   -value_type => 'ratio');

$dta->store($density_type);


my ( $current_start, $current_end, $current );
my $slice_count = 0;
my $block_size;

foreach my $slice (@sorted_slices){

  $block_size = $slice->length() / $bin_count;

  my @density_features=();

  print "GC percentage for ".$slice->seq_region_name().
    " with block size $block_size\n";

  $current_end = 0;
  $current = 0;

  while ($current_end < $slice->end()) {

    $current += $block_size;
    $current_start = $current_end+1;
    $current_end = int( $current + 1 );

    if ( $current_end < $current_start ) { 
      $current_end = $current_start;
    }

    if ( $current_end > $slice->end() ) {
      $current_end = $slice->end();
    }


    my $sub_slice = $slice->sub_Slice( $current_start , $current_end );

    my $gc = $sub_slice->get_base_count()->{'%gc'};
    my $df =  Bio::EnsEMBL::DensityFeature->new
      (-seq_region    => $slice,
       -start         => $current_start,
       -end           => $current_end,
       -density_type  => $density_type,
       -density_value => $gc);
    #print join ("\t", $slice, $current_start, $current_end, $density_type, $gc, "\n");
    $dfa->store($df);

  }

  last if ( $slice_count++ > $max_slices );
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




  


