use strict;

use Getopt::Long;
use XrefParser::BaseParser;

my ($host, $port, $dbname, $user, $pass, @species, @sources);

GetOptions('user=s'    => \$user,
	   'pass=s'    => \$pass,
	   'host=s'    => \$host,
	   'port=i'    => \$port,
	   'dbname=s'  => \$dbname,
	   'species=s' => \@species,
	   'source=s'  => \@sources,
	   'help'     => sub { usage(); exit(0); });

@species = split(/,/,join(',',@species));
@sources  = split(/,/,join(',',@sources));

if (!$user || !$host || !$dbname) {

  usage();
  exit(1);


}

XrefParser::BaseParser::run($host, $port, $dbname, $user, $pass, \@species, \@sources);

# --------------------------------------------------------------------------------

sub usage {

  print << "EOF";

  BaseParser.pm -user {user} -pass {password} -host {host} -port {port} -dbname {database} -species {species1,species2} -source {source1,source2}

EOF

}

# --------------------------------------------------------------------------------
