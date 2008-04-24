#!/software/bin/perl -w

#
# calculates the variation density from given core database
# It finds Variation database by itself using naming convention
# s/core/variation/
#
#
use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DensityType;
use Bio::EnsEMBL::DensityFeature;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Getopt::Long;

use Data::Dumper;
$Data::Dumper::Maxdepth = 2;

my $bin_count  = 150;
my $max_slices = 100;

my ( $host, $user, $pass, $port, $dbname  );

my ( $block_count, $genome_size, $block_size );

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
					    -dbname => $dbname,
					    -group => 'core',
					    -species => 'DEFAULT'
					   );

if( ! variation_attach( $db )) {
  die( "Couldnt attach variation to $dbname" );
}


#
# Clean up old features first. Also remove analysis and density type entry as these are recreated.
#

my $sth = $db->dbc->prepare("DELETE df, dt, a FROM density_feature df, density_type dt, analysis a WHERE a.analysis_id=dt.analysis_id AND dt.density_type_id=df.density_type_id AND a.logic_name='snpDensity'");
$sth->execute();

$sth = $db->dbc()->prepare(
  qq(
  DELETE ad
  FROM analysis_description ad
  WHERE ad.display_label = 'snpDensity') );
$sth->execute();

#
# Get the adaptors needed;
#

my $dfa = $db->get_DensityFeatureAdaptor();
my $dta = $db->get_DensityTypeAdaptor();
my $aa  = $db->get_AnalysisAdaptor();
my $slice_adaptor = $db->get_SliceAdaptor();


# Sort slices by coordinate system rank, then by length.
my @sorted_slices =
  sort( {      $a->coord_system()->rank() <=> $b->coord_system()->rank()
            || $b->seq_region_length() <=> $a->seq_region_length()
  } @{ $slice_adaptor->fetch_all('toplevel') } );

my $analysis =
  new Bio::EnsEMBL::Analysis(
              -program     => "variation_density.pl",
              -database    => "ensembl",
              -gff_source  => "variation_density.pl",
              -gff_feature => "density",
              -logic_name  => "snpDensity",
              -description => 'Density of SNP features on the sequence',
              -display_label => 'snpDensity',
              -displayable   => 1 );

$aa->store($analysis);
$aa->update($analysis);


#
# Create new density type.
#
my $dt = Bio::EnsEMBL::DensityType->new(-analysis   => $analysis,
					-region_features => $bin_count,
					-value_type => 'sum');
$dta->store($dt);


#
# Now the actual feature calculation loop
#

my ( $current_start, $current_end );

my $slice_count = 0;
my ( $current, $current_start, $current_end  );



foreach my $slice (@sorted_slices){

  $block_size = $slice->length() / $bin_count;


  print "SNP densities for ".$slice->seq_region_name().
    " with block size $block_size\n";

  $current_end = 0;
  $current = 0;

  while($current_end < $slice->end()) {

    $current += $block_size;
    $current_start = $current_end+1;
    $current_end = int( $current + 1 );

    if( $current_end < $current_start ) { 
      $current_end = $current_start;
    }

    if( $current_end > $slice->end() ) {
      $current_end = $slice->end();
    }


    my $sub_slice = $slice->sub_Slice( $current_start, $current_end );

    my $count =0;

    #
    #  How many snps fall into this subslice
    #
    foreach my $varf (@{$sub_slice->get_all_VariationFeatures()}){
      if( $varf->start >= 1 ) {
	$count++
      }
    }

    my $df = Bio::EnsEMBL::DensityFeature->new
      (-seq_region    => $slice,
       -start         => $current_start,
       -end           => $current_end,
       -density_type  => $dt,
       -density_value => $count);

    $dfa->store($df);
  }

  last if ( $slice_count++ > $max_slices );

}


#
# tries to attach variation database.
#

sub variation_attach {
  my $db = shift;

  my $core_db_name;
  $core_db_name = $db->dbc->dbname();
  if( $core_db_name !~ /_core_/ ) {
    return 0;
  }
  #
  # get a list of all databases on that server
  #
  my $sth = $db->dbc->prepare( "show databases" );
  $sth->execute();
  my $all_db_names = $sth->fetchall_arrayref();
  my %all_db_names = map {( $_->[0] , 1)} @$all_db_names;
  my $snp_db_name = $core_db_name;
  $snp_db_name =~ s/_core_/_variation_/;
 if( ! exists $all_db_names{ $snp_db_name } ) {
   return 0;
 }

 # this should register the dbadaptor with the Registry
 my $snp_db = Bio::EnsEMBL::Variation::DBSQL::DBAdaptor->new
   ( -host => $db->host(),
     -user => $db->username(),
     -pass => $db->password(),
     -port => $db->port(),
     -dbname => $snp_db_name,
     -group => "variation",
     -species => "DEFAULT"
   );

  return 1;
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
  if( !$max ) { $max = 1 };

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
