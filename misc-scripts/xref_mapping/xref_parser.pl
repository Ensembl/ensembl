use strict;

use Getopt::Long qw(:config pass_through);
use XrefParser::BaseParser;
use XrefParser::ProcessData;
use XrefParser::Database;
use XrefMapper::BasicMapper;

$| = 1;

my ( $host,             $port,          $dbname,
     $user,             $pass,          $species,
     $sources,          $checkdownload, $create,
     $release,          $cleanup,       $drop_existing_db,
     $deletedownloaded, $dl_path,       $notsource,
     $unzip,            $stats,         $notverbose,        $force, 
     $file,             $forceold );

my $options = join(" ",@ARGV);

print "Options: ".join(" ",@ARGV)."\n";

$unzip = 0;    # Do not decompress gzipped files by default

$notverbose = 0;
$port = 3306;

GetOptions(
    'file=s'         => \$file,
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
    'notverbose'     => \$notverbose,
    'notsource=s'    => \$notsource,
    'drop_db|dropdb!' => \$drop_existing_db,    # Drop xref db without user interaction
    'delete_downloaded' => \$deletedownloaded,
    'download_path=s' => \$dl_path,
    'unzip'           => \$unzip,                   # Force decompression of files
    'force'           => \$force,
    'forceold'        => \$forceold,
    'help'  => sub { usage(); exit(0); } );

if($ARGV[0]){
  print STDERR "Unknown command line arguments:-\n";
  foreach my $a (@ARGV){
    print STDERR "\t".$a."\n";
  }
  print STDERR "Stopping script. Please fix the command line.\n";
  print STDERR "use -help for full list of command line options.\n";;
  exit(1);
}

my @species;
my @sources  = split(/,/,join(',',$sources));
my @notsource  = split(/,/,join(',',$notsource));
my $dbc;
my $mapper;

if($file) {
  print "Using mapper configuration file\n";
  $mapper = XrefMapper::BasicMapper->process_file($file, !$notverbose);
  $dbc = $mapper->xref();
  #We do not support this in the new system
  if($species =~ /,/) {
    print STDERR "We do not support multiple species (-species human,mouse) and configuration files (-file) due to ambiguity in the target species required. Please do not use this option (it probably is not doing what you expect anyway)\n";
    exit 1;
  }
  elsif($species) {
    print STDERR "We do not pay any attention to -species when using the -file configuration files\n";
    exit 1;
  }
  @species = ($mapper->core()->species());
}
else {
  if(! $forceold) {
    print STDERR "Attempting to use non-config interface (some sources need a core DB so be warned). Please use a xref_mapper.input file or rerun with -forceold\n";
    exit 1;
  }
  print STDERR "Forced using old configuration; some sources e.g. RFAM need a core DB. Be warned\n";
  @species = split(/,/,join(',',$species));
  if ( !$user || !$host || !$dbname ) {
    usage();
    exit(1);
  }
  print "Host is $host\n";
  my $dbc = XrefParser::Database->new({
    host    => $host,
    dbname  => $dbname,
    port    => $port,
    user    => $user,
    pass    => $pass,
    verbose => (!$notverbose) 
  });
}

my $process  = XrefParser::ProcessData->new();

$process->run({ 
    dbc              => $dbc,
    mapper           => $mapper, # can go in undefined as we allow parsers to fail if they do not have access to the mapper
		speciesr         => \@species,
		sourcesr         => \@sources,
		checkdownload    => $checkdownload,
		create           => $create,
		release          => $release,
		cleanup          => $cleanup,
		drop_db          => $drop_existing_db,
		deletedownloaded => $deletedownloaded,
		dl_path          => (defined $dl_path ? $dl_path : "./"),
		notsourcesr      => \@notsource,
		unzip            => $unzip,
		stats            => $stats,
		verbose          => !($notverbose),
		force            => $force });

my $base_parser = XrefParser::BaseParser->new($process->database);

$base_parser->add_meta_pair("options",$options);


#
# If any sources are missing then we need to calculate the display xrefs from the core.
# As the xref database will not have all the data. Using the core is much slower!!
#

if($options =~ /source/ ){
  $base_parser->add_meta_pair("fullmode","no");
}
else{
  $base_parser->add_meta_pair("fullmode","yes");
}


$base_parser->parsing_finished_store_data();

# --------------------------------------------------------------------------------

sub usage {

  print << "EOF";

  xref_parser.pl -file {file} -user {user} -pass {password} -host {host} \\
    -port {port} -dbname {database} -species {species1,species2} \\
    -source {source1,source2} -notsource {source1,source2} \\
    -create -setrelease -deletedownloaded -checkdownload -stats -verbose \\
    -cleanup -drop_db -download_path -unzip \\
    -forceold

  -file             Use a xref_mapper.input configuration file to configure
                    this application. 

  -user             User name to access database. Must allow writing.
                    Ignored if using -file

  -pass             Password for user.
                    Ignored if using -file

  -host             Database host.
                    Ignored if using -file

  -port             Database port.
                    Ignored if using -file

  -dbname           Name of xref database to use/create.
                    Ignored if using -file

  -species          Which species to import. Multiple -species arguments
                    and/or comma, separated lists of species are
                    allowed. Species may be referred to by genus/species
                    (e.g. homo_sapiens) or common aliases (e.g. human).
                    Specifying an unknown species will cause a list
                    of valid species to be printed.  Not specifying a
                    -species argument will result in all species being
                    used.
                    
                    Ignored if using -file

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

  -notverbose       Do not output messages about the parsing. (NOT recomended)

  -stats            Generate the stats for the number of types of xrefs added.

  -force            No confirmation of actions just do it.
  
  -forceold			Force a run of the old syntax if required. This is fine so 
                long as you do not use a resource which expects to find a core
                database

EOF

}

#--------------------------------------------------------------------------------
