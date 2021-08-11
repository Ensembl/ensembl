#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2021] EMBL-European Bioinformatics Institute
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
use XrefParser::Database;
use Bio::EnsEMBL::Registry;
use File::Basename;

my ( $host,             $port,          $user,
     $pass,             $dbname,        $release,
     $species,          $taxon,         $division_id,
     $parser,           $source,        $file,
     $db,               $keep_db,       $help);

my $options = join(" ",@ARGV);

print "Options: ".join(" ",@ARGV)."\n";

GetOptions(
    'dbuser|user=s'  => \$user,
    'dbpass|pass=s'  => \$pass,
    'dbhost|host=s'  => \$host,
    'dbport|port=i'  => \$port,
    'dbname=s'       => \$dbname,
    'release=s'      => \$release,
    'species=s'      => \$species,
    'taxon=s'        => \$taxon,
    'division_id=s'  => \$division_id,
    'parser=s'       => \$parser,
    'source=s'       => \$source,
    'file=s'         => \$file,
    'db=s'           => \$db,
    'keep_db=s'      => \$keep_db,
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

  if ( !$host || !$species || !$parser) {
    usage();
    exit(1);
  }

  my $registry = 'Bio::EnsEMBL::Registry';
  $registry->load_registry_from_multiple_dbs(
    {
      -host       => $host,
      -port       => $port,
      -user       => $user,
      -pass       => $pass,
      -db_version => $release
    });

  if (!defined $taxon) {
    my $meta_container = $registry->get_adaptor($species,'core', 'MetaContainer');
    $taxon = $meta_container->get_taxonomy_id();
  }
  if (!defined $division_id) {
    my $meta_container = $registry->get_adaptor($species,'core', 'MetaContainer');
    my $division = $meta_container->get_division();
    my %division_taxon = (
    'Ensembl'            => 7742,
    'EnsemblVertebrates' => 7742,
    'Vertebrates'        => 7742,
    'EnsemblMetazoa'     => 33208,
    'Metazoa'            => 33208
    );
    $division_id = $division_taxon{$division};
  }

  my $sql_dir = dirname($0);

  my $xref_dbc = XrefParser::Database->new({
            host    => $host,
            dbname  => $dbname,
            port    => $port,
            user    => $user,
            pass    => $pass });
  $xref_dbc->create($sql_dir, 1, 1) unless $keep_db;
  my $xref_db_url = sprintf("mysql://%s:%s@%s:%s/%s", $user, $pass, $host, $port, $dbname);
  my $xref_dbi = $xref_dbc->dbi();

  my $module = "XrefParser::$parser";
  eval "require $module";
  my $xref_run = $module->new($xref_dbc);

  my $source_id = get_source_id($xref_dbi, $parser, $taxon, $source, $division_id);

  if (defined $db) {
    my $dba = $registry->get_DBAdaptor($species, $db);
    $dba->dbc()->disconnect_if_idle();
    $xref_run->run_script( { source_id  => $source_id,
                             species_id => $taxon,
                             dba        => $dba,
                             dbi        => $xref_dbi,
                             species    => $species,
                             file       => $file}) ;
  } else {
    my @files;
    push @files, $file;
    $xref_run->run( { source_id  => $source_id,
                      species_id => $taxon,
                      species    => $species,
                      dbi        => $xref_dbi,
                      files      => [@files] }) ;
  }

sub get_source_id {
  my ($dbi, $parser, $species_id, $name, $division_id) = @_;
  $name = "%$name%";
  my $source_id;
  my $select_source_id_sth = $dbi->prepare("SELECT u.source_id FROM source_url u, source s WHERE s.source_id = u.source_id AND parser = ? and species_id = ?");
  my $select_count_source_id_sth = $dbi->prepare("SELECT count(*) FROM source_url u, source s WHERE s.source_id = u.source_id AND parser = ? AND species_id = ?");
  $select_count_source_id_sth->execute($parser, $species_id);
  my $count = ($select_count_source_id_sth->fetchrow_array());
  if ($count == 1) {
    $select_source_id_sth->execute($parser, $species_id);
    $source_id = ($select_source_id_sth->fetchrow_array());
  }
  $select_source_id_sth = $dbi->prepare("SELECT u.source_id FROM source_url u, source s WHERE s.source_id = u.source_id AND parser = ? and species_id = ? and name like ?");
  $select_count_source_id_sth = $dbi->prepare("SELECT count(*) FROM source_url u, source s WHERE s.source_id = u.source_id AND parser = ? AND species_id = ? AND name like ?");
  $select_count_source_id_sth->execute($parser, $species_id, $name);
  $count = ($select_count_source_id_sth->fetchrow_array());
  if ($count == 1) {
    $select_source_id_sth->execute($parser, $species_id, $name);
    $source_id = ($select_source_id_sth->fetchrow_array());
  }
  # If no species-specific source, look for common sources
  if (!defined $source_id) {
    $select_source_id_sth->execute($parser, $division_id, $name);
    $source_id = ($select_source_id_sth->fetchrow_array())[0];
  }
  $select_source_id_sth->finish();
  $select_count_source_id_sth->finish();
  return $source_id;
}



# --------------------------------------------------------------------------------

sub usage {

  print << "EOF";

  xref_parser.pl -host {host} -port {port} -user {user} -pass {pass} -dbname {dbname} -release {release} \\
    -species {species} -taxon_id {taxon_id} \\
    -parser {parser} -source {source_id} -file {file} \\
    -db {db} =keep_db {keep_db} \\
    -help

  -user             User name to access database. Must allow writing.

  -pass             Password for user.

  -host             Database host.

  -port             Database port.

  -dbname           Name of xref database to use/create.

  -release          Release version of the species to parse
                    Used to find the right database on the server specified in the arguments

  -species          Which species to import. 
                    Species may be referred to by genus/species
                    (e.g. homo_sapiens) or common aliases (e.g. human).
                    Specifying an unknown species will cause a list
                    of valid species to be printed.

  -taxon            Which taxon to import.
                    Can be used as an alternative to species.

  -division         Which division the species belongs to.
                    This defines which sources will be parsed and does
                    not necessarily imply taxonomic relationship
                    (e.g. ciona intestinalis is a vertebrate in this context)

  -parser           Which parser to run

  -source           Name of the source to extra data for (should match equivalent parser)

  -file             Location and name of the file to be parsed
                    Path should be absolute

  -db               If the parser requires connection to a database, specify here
                    For example, specify otherfeatures when running RefSeqCoordinateParser

  -keep_db          When re-using an existing xref database, use the option
                    By default, deletes any existing one and creates a new one

EOF

}

#--------------------------------------------------------------------------------
