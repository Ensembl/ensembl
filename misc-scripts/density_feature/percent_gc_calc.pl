#!/usr/bin/perl -w

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


my $small_blocksize = 500;


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


#
# Get the adaptors needed;
#

my $slice_adaptor = $db->get_SliceAdaptor();
my $dfa = $db->get_DensityFeatureAdaptor();
my $dta = $db->get_DensityTypeAdaptor();
my $aa  = $db->get_AnalysisAdaptor();

my $slices = $slice_adaptor->fetch_all( "toplevel" );

my ( $large_blocksize, $genome_size );
for my $slice ( @$slices ) {
  $genome_size += $slice->length();
}

$large_blocksize = int( $genome_size / 4000 );


#
# Create new analysis object for density calculation.
#

my $analysis = new Bio::EnsEMBL::Analysis (-program     => "percent_gc_calc.pl",
					   -database    => "ensembl",
					   -gff_source  => "percent_gc_calc.pl",
					   -gff_feature => "density",
					   -logic_name  => "PercentGC");

$aa->store($analysis);


#
# Create new density type.
#

my $small_density_type = Bio::EnsEMBL::DensityType->new
  (-analysis   => $analysis,
   -block_size => $small_blocksize,
   -value_type => 'ratio');

my $large_density_type = Bio::EnsEMBL::DensityType->new
  (-analysis   => $analysis,
   -block_size => $large_blocksize,
   -value_type => 'ratio');

$dta->store($small_density_type);
$dta->store($large_density_type);


my ( $current_start, $current_end );

foreach my $slice ( @$slices ) {

#
# do it for small and large blocks
#

  for my $density_type ( $large_density_type, $small_density_type ) {

    my $blocksize = $density_type->block_size();
    $current_start = 1;

    while($current_start <= $slice->end()) {
      $current_end = $current_start+$blocksize-1;
      if( $current_end > $slice->end() ) {
        $current_end = $slice->end();
      }

      my $sub_slice = $slice->sub_Slice( $current_start, $current_end );

      my $gc = $sub_slice->get_base_count()->{'%gc'};
      my $df =  Bio::EnsEMBL::DensityFeature->new
        (-seq_region    => $slice,
         -start         => $current_start,
         -end           => $current_end,
         -density_type  => $density_type,
         -density_value => $gc);

      $dfa->store($df);

      $current_start = $current_end+1;
    }
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




  


