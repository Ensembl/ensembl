#
# script to calculate the snp density with help of
#  attached lite database. Only works if argument database
#  is a core database and lite can be found by substituting
#  s/_core_/_lite_/
#
# blocksize condition is 4_000 per genome?


use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DensityType;
use Bio::EnsEMBL::DensityFeature;
use Bio::EnsEMBL::Lite::DBAdaptor;
use Getopt::Long;

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
					    -dbname => $dbname);

if( ! lite_attach( $db )) {
  die( "Couldnt attache lite to $dbname" );
}

#
# Get the adaptors needed;
#

my $dfa = $db->get_DensityFeatureAdaptor();
my $dta = $db->get_DensityTypeAdaptor();
my $aa  = $db->get_AnalysisAdaptor();
my $slice_adaptor = $db->get_SliceAdaptor();

my $top_slices = $slice_adaptor->fetch_all( "toplevel" );


my ( $block_size, $genome_size );
for my $slice ( @$top_slices ) {
  $genome_size += $slice->length();
}

$block_size = int( $genome_size / 4000 );

my $analysis = new Bio::EnsEMBL::Analysis (-program     => "snp_lite_density.pl",
					   -database    => "ensembl",
					   -gff_source  => "snp_lite_density.pl",
					   -gff_feature => "density",
					   -logic_name  => "snpDensity");

$aa->store( $analysis );


#
# Create new density type.
#
my $dt = Bio::EnsEMBL::DensityType->new(-analysis   => $analysis,
					-block_size => $block_size,
					-value_type => 'sum');
$dta->store($dt);


#
# Now the actual feature calculation loop
#

my ( $current_start, $current_end );



foreach my $slice (@$top_slices){

  $current_start = 1;

  print "SNP densities for ".$slice->seq_region_name().
    " with block size $block_size\n";

  while($current_start <= $slice->end()) {
    $current_end = $current_start+$block_size-1;
    if( $current_end > $slice->end() ) {
      $current_end = $slice->end();
    }


    my $sub_slice = $slice->sub_Slice( $current_start, $current_end );

    my $count =0;

    #
    #  How many snps fall into this subslice
    #
    foreach my $snp (@{$sub_slice->get_all_SNPs()}){
      if( $snp->start >= 1 ) {
        $count++
      }
    }

    my $df = Bio::EnsEMBL::DensityFeature->new
      (-seq_region    => $slice,
       -start         => $current_start,
       -end           => $current_end,
       -density_type  => $dt,
       -density_value => $count);

    $current_start = $current_end + 1;

    $dfa->store($df);
  }
}


#
# tries to attach lite.
#

sub lite_attach {
  my $db = shift;

  my $core_db_name;
  $core_db_name = $db->dbname();
  if( $core_db_name !~ /_core_/ ) {
    return 0;
  }
  #
  # get a lost of all databases on that server
  #
  my $sth = $db->prepare( "show databases" );
  $sth->execute();
  my $all_db_names = $sth->fetchall_arrayref();
  my %all_db_names = map {( $_->[0] , 1)} @$all_db_names;
  my $snp_db_name = $core_db_name;
  $snp_db_name =~ s/_core_/_lite_/;
  if( ! exists $all_db_names{ $snp_db_name } ) {
    return 0;
  }

  my $snp_db = Bio::EnsEMBL::Lite::DBAdaptor->new
    ( -host => $db->host(),
      -user => $db->username(),
      -pass => $db->password(),
      -port => $db->port(),
      -dbname => $snp_db_name );
  $db->add_db_adaptor( "lite", $snp_db );
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




  


