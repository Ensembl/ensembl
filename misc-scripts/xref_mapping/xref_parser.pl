#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

use strict;
use warnings;

use Getopt::Long qw(:config pass_through);
use XrefParser::BaseParser;
use XrefParser::ProcessData;

my ( $host,             $port,          $dbname,
     $user,             $pass,          $species,
     $sources,          $checkdownload, $create,
     $release,          $cleanup,       $drop_existing_db,
     $deletedownloaded, $dl_path,       $notsource,
     $unzip, $stats, $notverbose, $force );

my $options = join(" ",@ARGV);

print "Options: ".join(" ",@ARGV)."\n";

$unzip = 0;    # Do not decompress gzipped files by default

$notverbose = 0;

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
    'notverbose'     => \$notverbose,
    'notsource=s'    => \$notsource,
    'drop_db|dropdb!' => \$drop_existing_db,    # Drop xref db without user interaction
    'delete_downloaded' => \$deletedownloaded,
    'download_path=s' => \$dl_path,
    'unzip'           => \$unzip,                   # Force decompression of files
    'force'           => \$force,
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

my @species = split(/,/,join(',',$species)) if $species;
my @sources  = split(/,/,join(',',$sources)) if $sources;
my @notsource  = split(/,/,join(',',$notsource)) if $notsource;

$| = 1;

if ( !$user || !$host || !$dbname ) {
    usage();
    exit(1);
}

if ( $host =~ /staging/ || $host =~ /livemirror/ ){
  print STDERR "$host was specified as host name. This is a production server and should not be used to create the xref database\n";
  exit(1);
}

if ($dbname =~ /_core_/) {
  print STDERR "$dbname was specified for the database name. This should be the name of the xref database, not the core database\n";
  exit(1);
}


print "host os $host\n";

my $process  = XrefParser::ProcessData->new();

$process->run({ host             => $host,
		port             =>  ( defined $port ? $port : '3306' ),
		dbname           => $dbname,
		user             => $user,
		pass             => $pass,
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

  -notverbose       Do not output messages about the parsing. (NOT recomended)

  -stats            Generate the stats for the number of types of xrefs added.

  -force            No confirmation of actions just do it.

EOF

}

#--------------------------------------------------------------------------------
