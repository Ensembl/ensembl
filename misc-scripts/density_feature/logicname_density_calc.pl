#
# script to calculate the density of features 
# for any given analysis.logic_name
#

use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DensityType;
use Bio::EnsEMBL::DensityFeature;
use Getopt::Long;

my ( $host, $user, $pass, $port, $dbname, $logic_name );

my ( $block_count, $genome_size, $block_size );

GetOptions( "host=s", \$host,
	    "user=s", \$user,
	    "pass=s", \$pass,
	    "port=i", \$port,
	    "dbname=s", \$dbname,
            "logic_name=s",\$logic_name,
	  );

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $host,
					    -user => $user,
					    -port => $port,
					    -pass => $pass,
					    -dbname => $dbname);

$logic_name || die( "Need --logic_name arg" );

my $analysis_adaptor = $db->get_adaptor("Analysis");
my $analysis = $analysis_adaptor->fetch_by_logic_name( $logic_name ) || 
    die( "logic_name $logic_name not found in database" );

#
# What sort of feature does this logic_name correspond to?
# Need this to estimate block size etc
#
my @feature_types;
foreach my $type( $analysis_adaptor->feature_classes() ){
  foreach my $analysis( @{$analysis_adaptor->fetch_all_by_feature_class($type)} ){
    if( uc($analysis->logic_name) eq uc($logic_name) ){
      push( @feature_types, $type );
      last;
    }
  }
}
if( ! @feature_types ){
  die( "Logic name $logic_name is a valid analysis in the database, ".
       "but does not correspond to any features\n" );
}
if( @feature_types > 1 ){
  # Hmm - logic name maps to more than one feature type!
  warn( "$logic_name maps to multiple feature types: ".
        join(', ', @feature_types ) );
}
my $feature_type = $feature_types[0];
warn( "Logic_name $logic_name maps to features of type $feature_type\n" );

#
# Check that database holds an assembly
#
my $sth = $db->dbc()->prepare( "select count(*)  from seq_region" );
$sth->execute();
my ( $seq_region_count ) = $sth->fetchrow_array();
$seq_region_count > 0 or die( "No seq_regions for $dbname\n" );


#
# Get the adaptors needed;
#

my $dfa           = $db->get_adaptor('DensityFeature');
my $dta           = $db->get_adaptor('DensityType');
my $slice_adaptor = $db->get_adaptor('Slice');

#
# block size estimation
#
my $feature_table = join( '_', map lc, ( $feature_type =~ /([A-Z][a-z]+)/g ) );
my $analysis_id = $analysis->dbID;
my $count_sql = qq(
SELECT COUNT(*) 
FROM  $feature_table 
WHERE analysis_id=$analysis_id );
my $sth = $db->dbc()->prepare( $count_sql );
$sth->execute() || die( $sth->errstr );
my $feature_count = $sth->fetchrow_array();
warn( "$feature_count features found for $logic_name\n" );
$block_count = $feature_count >> 1;

my $top_slices = $slice_adaptor->fetch_all('toplevel');
for my $slice ( @$top_slices ) {
  $genome_size += $slice->length();
}

$block_size = int( $genome_size / $block_count );
	

#
# Create new analysis entry based lon logic_name
#
$analysis = new Bio::EnsEMBL::Analysis 
    (-program     => "logic_name_calc.pl",
     -database    => "ensembl",
     -gff_source  => "logic_name_calc.pl",
     -gff_feature => "density",
     -logic_name  => "density_$logic_name");
$analysis_adaptor->store( $analysis );

#
# Now the actual feature calculation
#
# Create new density type.
#
$analysis = $analysis_adaptor->fetch_by_logic_name("density_$logic_name");
my $dt = Bio::EnsEMBL::DensityType->new(-analysis   => $analysis,
                                        -block_size => $block_size,
                                        -value_type => 'sum');
$dta->store($dt);

my ( $current_start, $current_end );

my $class = join('', map{ucfirst($_)} split('_', $feature_type) );
my $getter_method = "get_all_${class}s";

foreach my $slice (@$top_slices){

  $current_start = 1;

  my @density_features=();

  print "$logic_name densities for ".$slice->seq_region_name().
      " with block size $block_size\n";

  while($current_start <= $slice->end()) {
    $current_end = $current_start+$block_size-1;
    if( $current_end > $slice->end() ) {
      $current_end = $slice->end();
    }
    
    my $sub_slice = $slice->sub_Slice( $current_start, $current_end );
    
    my $count = scalar( @{$sub_slice->$getter_method($logic_name)} );
    push @density_features, Bio::EnsEMBL::DensityFeature->new
        (-seq_region    => $slice,
	 -start         => $current_start,
	 -end           => $current_end,
	 -density_type  => $dt,
	 -density_value => $count);
    
    $current_start = $current_end + 1;
  }
  $dfa->store(@density_features);
  print "Created ", scalar @density_features, " $logic_name density features.\n";
  # print_features(\@density_features);
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




  


