use strict;

use Getopt::Long;
use XrefParser::BaseParser;

my ( $host,             $port,          $dbname,
     $user,             $pass,          $species,
     $sources,          $checkdownload, $create,
     $release,          $cleanup,       $drop_existing_db,
     $deletedownloaded, $dl_path,       $notsource,
     $unzip, $stats, $verbose );

print "Options: ".join(" ",@ARGV)."\n";

$unzip = 0;    # Do not decompress gzipped files by default

GetOptions(
    'dbuser|user=s'  => \$user,
    'dbpass|pass=s'  => \$pass,
    'dbhost|host=s'  => \$host,
    'dbport|port=i'  => \$port,
    'dbname=s'       => \$dbname,
    'species=s'      => \$species,
    'source=s'       => \$sources,
    'download_dir=s' => \$dl_path,
    'checkdownload!' => \$checkdownload,    # Don't download if exists
    'create'         => \$create,
    'setrelease=s'   => \$release,
    'cleanup'        => \$cleanup,
    'stats'          => \$stats,
    'verbose'        => \$verbose,
    'notsource=s'    => \$notsource,
    'drop_db|dropdb!' =>
      \$drop_existing_db,    # Drop xref db without user interaction
    'delete_downloaded' => \$deletedownloaded,
    'download_path=s'   => \$dl_path,
    'unzip' => \$unzip,                   # Force decompression of files
    'help'  => sub { usage(); exit(0); } );

my @species = split(/,/,join(',',$species));
my @sources  = split(/,/,join(',',$sources));
my @notsource  = split(/,/,join(',',$notsource));

$| = 1;

if ( !$user || !$host || !$dbname ) {
    usage();
    exit(1);
}

my $base_parser = XrefParser::BaseParser->new();

$base_parser->run(
               $host, ( defined $port ? $port : '3306' ),
               $dbname,           $user,
               $pass,             \@species,
               \@sources,         $checkdownload,
               $create,           $release,
               $cleanup,          $drop_existing_db,
               $deletedownloaded, $dl_path,
               \@notsource,       $unzip, $stats, $verbose );

# --------------------------------------------------------------------------------

sub usage {

  print << "EOF";

  xref_parser.pl -user {user} -pass {password} -host {host} \\
    -port {port} -dbname {database} -species {species1,species2} \\
    -source {source1,source2} -notsource {source1,source2} \\
    -create -setrelease -deletedownloaded -checkdownload -stats -verbose \\
    -cleanup -drop_db -download_path -unzip

  -user             User name to access database. Must allow writing.

  -pass             Password for user.

  -host             Database host.

  -port             Database port.

  -dbname           Name of xref database to use/create.

  -species          Which species to import. Multiple -species arguments
                    and/or comma, separated lists of species are
                    allowed. Species may be referred to by genus/species
                    (e.g. homo_sapiens) or common aliases (e.g. human).
                    Specifying an unknown species will cause a list
                    of valid species to be printed.  Not specifying a
                    -species argument will result in all species being
                    used.

  -source           Which sources to import. Multiple -source arguments
                    and/or comma, separated lists of sources are
                    allowed.  Specifying an unknown source will cause a
                    list of valid sources to be printed.  Not specifying
                    a -source argument will result in all species being
                    used.

  -notsource        Which source to skip.

  -create           If specified, cause dbname to be deleted and
                    re-created if it already exists.  The user is
                    prompted before the database is dropped to
                    prevent disasters arising from dropping the wrong
                    database.  If needed, this switch will also make
                    the script ask the user whether it should run the
                    'xref_config2sql.pl' script for him/her to renew
                    'sql/populate_metadata.sql' from 'xref_config.ini'.

  -checkdownload    Check if file exists, otherwise downloads the file.

  -deletedownloaded Delete any existing downloaded files first.

  -setrelease       Set the release version for ALL the sources specified.

  -cleanup          Delete the downloaded source files after parsing.

  -drop_db          Drop the xref-database without user interaction.

  -download_path    Directory into which to download files (default is
                    the current directory).

  -unzip            Decompress gzipped files (default is to use compressed
                    files).

  -verbose          Output messages about the parsing.

  -stats            Generate the stats for the number of types of xrefs added.

EOF

}

#--------------------------------------------------------------------------------
