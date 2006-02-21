# Load data from a file of regulatory regions into a database

use strict;

use DBI;
use RegulatoryFeatureParser::BaseParser;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;

my ($host, $user, $pass, $port, $dbname, $file, $del, $type);

GetOptions( "host=s",   \$host,
	    "user=s",   \$user,
	    "pass=s",   \$pass,
	    "port=i",   \$port,
	    "dbname=s", \$dbname,
	    "type=s",   \$type,
	    "file=s",   \$file,
	    "delete",   \$del,
	    "help",     \&usage
	  );

my $db_adaptor = new Bio::EnsEMBL::DBSQL::DBAdaptor('-host' => $host,
						    '-port' => $port,
						    '-user' => $user,
						    '-pass' => $pass,
						    '-dbname' => $dbname);

usage() if (!$host || !$user || !$dbname || !$type);

# validate type
exit(1) if (!RegulatoryFeatureParser::BaseParser::validate_type($db_adaptor, $type));

RegulatoryFeatureParser::BaseParser::delete_existing($db_adaptor, $type) if ($del);


# load module corresponding to type
my $module;
eval "require RegulatoryFeatureParser::$type";
if($@) {
  warn("Did not find a parser module corresponding to $type, exiting.\n");
  exit(1);
} else{
  $module = $type;
}

my $parser = "RegulatoryFeatureParser::$type"->new();

if ($file) {

  die "Can't find $file\n" if (!-e $file);

  my $objects = $parser->parse($db_adaptor, $file);

  $parser->upload_features_and_factors($db_adaptor, $objects);

}




sub usage {
  print <<EOF

Usage: perl $0 <options>

  -host    Database host to connect to

  -port    Database port to connect to

  -dbname  Database name to use

  -user    Database username

  -pass    Password for user

  -type    Type of regulatory features to upload; must be present in analysis table and
           have a correspoinding parser

  -file    File to load data from; must be parseable by parser corresponding to type

  -delete  Delete all data related to regulatory features of specified type before loading

  -help    This message

EOF

}

