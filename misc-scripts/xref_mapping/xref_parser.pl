use strict;

use Getopt::Long;
use XrefParser::BaseParser;

my ($host, $port, $dbname, $user, $pass, @species, @sources, $skipdownload, $create, $release);

GetOptions('user=s'       => \$user,
	   'pass=s'       => \$pass,
	   'host=s'       => \$host,
	   'port=i'       => \$port,
	   'dbname=s'     => \$dbname,
	   'species=s'    => \@species,
	   'source=s'     => \@sources,
	   'skipdownload' => \$skipdownload,
	   'create'       => \$create,
	   'setrelease=s' => \$release,
	   'help'         => sub { usage(); exit(0); });

@species = split(/,/,join(',',@species));
@sources  = split(/,/,join(',',@sources));


if (!$user || !$host || !$dbname) {

  usage();
  exit(1);

}

XrefParser::BaseParser::run($host, $port, $dbname, $user, $pass, \@species, \@sources, $skipdownload, $create, $release);

# --------------------------------------------------------------------------------

sub usage {

  print << "EOF";

  xref_parser.pl -user {user} -pass {password} -host {host} -port {port} -dbname {database} -species {species1,species2} -source {source1,source2} -skipdownload -create -setrelease

  -user         User name to access database. Must allow writing.

  -pass         Password for user.

  -host         Database host.

  -port         Database port.

  -dbname       Name of xref database to use/create.

  -species      Which species to import. Multiple -species arguments and/or comma,
                separated lists of species are allowed. Species may be referred to
                by genus/species (e.g. homo_sapiens) or common aliases (e.g. human).
                Specifying an unknown species will cause a list of valid species to
                be printed.
                Not specifying a -species argument will result in all species being
                used.

  -source       Which sources to import. Multiple -source arguments and/or comma,
                separated lists of sources are allowed.
                Specifying an unknown source will cause a list of valid sources to
                be printed.
                Not specifying a -source argument will result in all species being
                used.

  -create       If specified, cause dbname to be deleted and re-created if it
                already exists.

  -skipdownload Don't download new data, parse existing.

  -setrelease   Set the release version for ALL the sources specified.

EOF

}

#--------------------------------------------------------------------------------
