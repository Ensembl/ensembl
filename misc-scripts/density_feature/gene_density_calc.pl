#
# script to calculate the gene density features on a database
# should work on any species database
#

#
# It will only run on databases with genes ...
# boundary condition: on average there should be 2 genes per block
#

use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DensityType;
use Bio::EnsEMBL::DensityFeature;
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

my $sth = $db->dbc()->prepare( "select count(*) from gene" );
$sth->execute();

my ( $gene_count )  = $sth->fetchrow_array();

if( ! $gene_count ) {
  print STDERR "No gene density for $dbname.\n";
  exit();
} else {
  $block_count = $gene_count >> 1;
}

#
# Could be database without seq_regions
#  Then have to try and attach core db
#
$sth = $db->dbc()->prepare( "select count(*)  from seq_region" );
$sth->execute();
my ( $seq_region_count ) = $sth->fetchrow_array();
if( ! $seq_region_count ) {
  #
  # for the time being only do core dbs
  # no dbs with no seq regions
  #
  print STDERR "No gene density for $dbname, no seq_regions.\n";
  exit();

  if( ($dbname =~ /_estgene_/ ) || ( $dbname =~ /_vega_/ )) {
    my $dna_db_name = $dbname;
    $dna_db_name =~ s/(_estgene_|_vega_)/_core_/;
    my $dna_db =  new Bio::EnsEMBL::DBSQL::DBAdaptor
      (-host => $host,
       -user => $user,
       -port => $port,
       -pass => $pass,
       -dbname => $dna_db_name );
    print STDERR "Attaching $dna_db_name to $dbname.\n";
    $db->dnadb( $dna_db );
  } else {
    print STDERR "No gene density for $dbname, no seq_regions.\n";
    exit();
  }
}

#
# Get the adaptors needed;
#

my $dfa = $db->get_DensityFeatureAdaptor();
my $dta = $db->get_DensityTypeAdaptor();
my $aa  = $db->get_AnalysisAdaptor();
my $slice_adaptor = $db->get_SliceAdaptor();

#
# block size estimation
#

my $top_slices = $slice_adaptor->fetch_all('toplevel');
for my $slice ( @$top_slices ) {
  $genome_size += $slice->length();
}

$block_size = int( $genome_size / $block_count );
	
my $analysis = new Bio::EnsEMBL::Analysis (-program     => "gene_density_calc.pl",
					   -database    => "ensembl",
					   -gff_source  => "gene_density_calc.pl",
					   -gff_feature => "density",
					   -logic_name  => "knownGeneDensity");

$aa->store( $analysis );

$analysis = new Bio::EnsEMBL::Analysis (-program     => "gene_density_calc.pl",
					-database    => "ensembl",
					-gff_source  => "gene_density_calc.pl",
					-gff_feature => "density",
					-logic_name  => "geneDensity");

$aa->store( $analysis );

#
# Now the actual feature calculation loop
#


foreach my $known (1, 0) {
  #
  # Create new analysis object for density calculation.
  #

  if($known) {
    $analysis = $aa->fetch_by_logic_name('knownGeneDensity');
  } else {
    $analysis = $aa->fetch_by_logic_name('geneDensity');
  }

  #
  # Create new density type.
  #

  my $dt = Bio::EnsEMBL::DensityType->new(-analysis   => $analysis,
					  -block_size => $block_size,
					  -value_type => 'sum');

  $dta->store($dt);

  my ( $current_start, $current_end );

  foreach my $slice (@$top_slices){

    $current_start = 1;

    my @density_features=();

    print "Gene densities for ".$slice->seq_region_name().
      " with block size $block_size\n";

    while($current_start <= $slice->end()) {
      $current_end = $current_start+$block_size-1;
      if( $current_end > $slice->end() ) {
	$current_end = $slice->end();
      }


      my $sub_slice = $slice->sub_Slice( $current_start, $current_end );

      my $count =0;

      #
      # Store info for genes (ignore pseudo genes)
      #

      foreach my $gene (@{$sub_slice->get_all_Genes()}){
	if($gene->type() !~ /pseudogene/i and $gene->start >=1 ) {
	  $count++ if(!$known || $gene->is_known());
	}
      }

      push @density_features, Bio::EnsEMBL::DensityFeature->new
	(-seq_region    => $slice,
	 -start         => $current_start,
	 -end           => $current_end,
	 -density_type  => $dt,
	 -density_value => $count);

      $current_start = $current_end + 1;
#      print STDERR ".";
    }
    $dfa->store(@density_features);
    print "Created ", scalar @density_features, " gene density features.\n";
    # print_features(\@density_features);
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




  


