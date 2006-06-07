use strict;

use Getopt::Long;
use XrefParser::BaseParser;

my ($host, $port, $dbname, $user, $pass, @species, @sources, $skipdownload, $checkdownload, $create, $release, $cleanup, $drop_existing_db, $deletedownloaded, $dl_path);

GetOptions('dbuser|user=s'       => \$user,
	   'dbpass|pass=s'       => \$pass,
	   'dbhost|host=s'       => \$host,
	   'dbport|port=i'       => \$port,
	   'dbname=s'            => \$dbname,
	   'species=s'           => \@species,
	   'source=s'            => \@sources,
	   'download_dir=s'      => \$dl_path,
	   'delete_downloaded'   => \$delete_downloaded ,  # deletes all downloaded files  
	   'skipdownload'        => \$skipdownload,     # skips all downloads
	   'checkdownload!'      => \$checkdownload,   # if file exists it won't be downloaded 
	   'create'       => \$create,
	   'setrelease=s' => \$release,
	   'cleanup'      => \$cleanup,
	   'drop_db|dropdb!'     => \$drop_existing_db, # drops xref db without user interaction
	   'delete_downloaded' => \$deletedownloaded,
	   'download_path=s' => \$dl_path,
	   'help'         => sub { usage(); exit(0); });

@species = split(/,/,join(',',@species));
@sources  = split(/,/,join(',',@sources));


if (!$user || !$host || !$dbname) {

  usage();
  exit(1);

}

XrefParser::BaseParser::run($host, $port, $dbname, $user, $pass, \@species, \@sources, $skipdownload, $checkdownload, $create, $release, $cleanup,$drop_existing_db, $deletedownloaded, $dl_path);

# --------------------------------------------------------------------------------

sub usage {

  print << "EOF";

  xref_parser.pl -user {user} -pass {password} -host {host} -port {port} -dbname {database} -species {species1,species2} -source {source1,source2} -skipdownload -create -setrelease

  -user             User name to access database. Must allow writing.
		
  -pass             Password for user.
		
  -host             Database host.
		
  -port             Database port.
		
  -dbname           Name of xref database to use/create.
		
  -species          Which species to import. Multiple -species arguments and/or comma,
                    separated lists of species are allowed. Species may be referred to
                    by genus/species (e.g. homo_sapiens) or common aliases (e.g. human).
                    Specifying an unknown species will cause a list of valid species to
                    be printed.
                    Not specifying a -species argument will result in all species being
                    used.
		
  -source           Which sources to import. Multiple -source arguments and/or comma,
                    separated lists of sources are allowed.
                    Specifying an unknown source will cause a list of valid sources to
                    be printed.
                    Not specifying a -source argument will result in all species being
                    used.
		
  -create           If specified, cause dbname to be deleted and re-created if it
                    already exists. User is prompted before database is dropped to
                    prevent disasters arising from dropping the wrong database.
		
  -skipdownload     Don't download any data, parse existing.

  -checkdownload    Check if file exists, otherwise downloads the file 

  -deletedownloaded Delete any existing downloaded files first.

  -setrelease       Set the release version for ALL the sources specified.

  -cleanup          Delete the downloaded source files after parsing.

  -drop_db         Drop the xref-database without user interaction

  -download_path   Directory into which to download files. Default is current directory.

EOF

}

#--------------------------------------------------------------------------------
