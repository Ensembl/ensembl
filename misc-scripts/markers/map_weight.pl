# script to calculate map_weights in a database that has markers
# and marker_features. Recreates the marker_feature table with weights set

use strict;
use DBI;

use Getopt::Long;

my ( $host, $user, $pass, $port, $dbname );
my $verbose = 0;

GetOptions( "host=s", \$host,
	    "user=s", \$user,
	    "pass=s", \$pass,
	    "port=i", \$port,
	    "dbname=s", \$dbname,
	    "verbose", \$verbose
	  );

if( !$host ) {
  usage();
}



my $dsn = "DBI:mysql:host=$host;dbname=$dbname";
if( $port ) {
  $dsn .= ";port=$port";
}

my $db = DBI->connect( $dsn, $user, $pass );


$db->do( "
  CREATE TABLE tmp_m_weight
  SELECT marker_id, count(*) as count 
  FROM marker_feature
  GROUP BY marker_id
" );

$db->do( "
  CREATE TABLE new_marker_feature
  SELECT mf.marker_feature_id, mf.marker_id, mf.seq_region_id, mf.seq_region_start,
         mf.seq_region_end, mf.analysis_id, tmw.count
  FROM   marker_feature mf, tmp_m_weight tmw
  WHERE  mf.marker_id = tmw.marker_id
" );

$db->do( "delete from marker_feature" );
$db->do( "insert into marker_feature select * from new_marker_feature" );
$db->do( "drop table tmp_m_weight" );
$db->do( "drop table new_marker_feature" );

sub usage {
  print <<EOF;
    
Usage: perl map_weight.pl [options]
   -user username for a write enabled user
   -host hostname
   -port portnumber
   -pass password
   -dbname database name where the markers and the features are

EOF

  exit;
}
