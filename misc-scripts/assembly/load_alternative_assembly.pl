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


=head1 NAME

load_alternative_assembly.pl - create a db for transfering annotation to the
Ensembl assembly

=head1 SYNOPSIS

load_alternative_assembly.pl [arguments]

Required arguments:

    --dbname, db_name=NAME              database name NAME
    --host, --dbhost, --db_host=HOST    database host HOST
    --port, --dbport, --db_port=PORT    database port PORT
    --user, --dbuser, --db_user=USER    database username USER
    --pass, --dbpass, --db_pass=PASS    database passwort PASS
    --assembly=ASSEMBLY                 assembly version ASSEMBLY

    --altdbname=NAME                    alternative database NAME
    --altassembly=ASSEMBLY              alternative assembly version ASSEMBLY
    --coord_systems, --cs=CS            coord_systems to load (default: all
                                        toplevel)

Optional arguments:

    --conffile, --conf=FILE             read parameters from FILE
                                        (default: conf/Conversion.ini)

    --logfile, --log=FILE               log to FILE (default: *STDOUT)
    --logpath=PATH                      write logfile to PATH (default: .)
    --logappend, --log_append           append to logfile (default: truncate)

    -v, --verbose=0|1                   verbose logging (default: false)
    -i, --interactive=0|1               run script interactively (default: true)
    -n, --dry_run, --dry=0|1            dont write results to database
    -h, --help, -?                      print help (this message)

=head1 DESCRIPTION

This script is part of a series of scripts to create a mapping between two
assemblies. It assembles the toplevel coordinate systems of two different
assemblies of a genome by creating a whole genome alignment between the two.

The process assumes that the two assemblies are reasonably similar, i.e. there
are no major rearrangements or clones moved from one toplevel seq_region to
another.

See "Related files" below for an overview of the whole process.

This particular script loads the alternative toplevel seq_regions into the
Ensembl database for further processing.

=head1 RELATED FILES

The whole process of creating a whole genome alignment between two assemblies
is done by a series of scripts. Please see

  ensembl/misc-scripts/assembly/README

for a high-level description of this process, and POD in the individual scripts
for the details.


=head1 AUTHOR

Patrick Meidl <meidl@ebi.ac.uk>, Ensembl core API team

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<http://lists.ensembl.org/mailman/listinfo/dev>

=cut

use strict;
use warnings;
no warnings 'uninitialized';

use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Utils::ConversionSupport;
use Cwd;


$| = 1;

my $path = getcwd;

my $support = new Bio::EnsEMBL::Utils::ConversionSupport($path);

# parse options
$support->parse_common_options(@_);
$support->parse_extra_options(
    'assembly=s',
    'altdbname=s',
    'altassembly=s',
    'coord_systems|cs=s',
);
$support->allowed_params(
    $support->get_common_params,
    'assembly',
    'altdbname',
    'altassembly',
    'coord_systems',
);

if ($support->param('help') or $support->error) {
    warn $support->error if $support->error;
    pod2usage(1);
}

$support->comma_to_list('coord_systems');

# ask user to confirm parameters to proceed
$support->confirm_params;

# get log filehandle and print heading and parameters to logfile
$support->init_log;

$support->check_required_params(
  'assembly',
  'altdbname',
  'altassembly'
);

if ($support->param('dry_run')) {
  $support->log("Nothing to do for a dry run. Exiting.\n\n");
  $support->finish_log;
  exit;
}

#####
# connect to database and get adaptors
#
my ($dba, $dbh, $sql, $sth);

# first set connection parameters for alternative db
# both databases have to be on the same host, so we don't need to configure
# them separately
map { $support->param("alt$_", $support->param($_)) } qw(host port user pass);

# reference database
$dba->{'ref'} = $support->get_database('ensembl');
$dbh->{'ref'} = $dba->{'ref'}->dbc->db_handle;

# database containing the alternative assembly
$dba->{'alt'} = $support->get_database('core', 'alt');
$dbh->{'alt'} = $dba->{'alt'}->dbc->db_handle;

#####
# create backups of the tables that will be modified
#
$support->log_stamped("Creating table backups...\n");
$support->log_stamped("seq_region...\n", 1);
$dbh->{'ref'}->do('CREATE TABLE seq_region_bak LIKE seq_region');
$dbh->{'ref'}->do('INSERT INTO seq_region_bak SELECT * FROM seq_region');
$support->log_stamped("seq_region_attrib...\n", 1);
$dbh->{'ref'}->do('CREATE TABLE seq_region_attrib_bak LIKE seq_region_attrib');
$dbh->{'ref'}->do('INSERT INTO seq_region_attrib_bak SELECT * FROM seq_region_attrib');
$support->log_stamped("assembly...\n", 1);
$dbh->{'ref'}->do('CREATE TABLE assembly_bak LIKE assembly');
$dbh->{'ref'}->do('INSERT INTO assembly_bak SELECT * FROM assembly');
$support->log_stamped("meta...\n", 1);
$dbh->{'ref'}->do('CREATE TABLE meta_bak LIKE meta');
$dbh->{'ref'}->do('INSERT INTO meta_bak SELECT * FROM meta');
$support->log_stamped("coord_system...\n", 1);
$dbh->{'ref'}->do('CREATE TABLE coord_system_bak LIKE coord_system');
$dbh->{'ref'}->do('INSERT INTO coord_system_bak SELECT * FROM coord_system');
$support->log_stamped("Done.\n\n");

#####
# load seq_regions from alternative assembly db
#
$support->log_stamped("Load seq_regions from alternative db...\n");

# determine which coord_systems we want to include
# by default, all coord_systems with toplevel seq_regions will be used
my @coord_systems = $support->param('coord_systems');

unless (@coord_systems) {
  # get toplevel coord_systems from both dbs
  my %common_cs = ();

  my $sql = qq(
    SELECT distinct(cs.name)
    FROM seq_region sr, coord_system cs, seq_region_attrib sra, attrib_type at
    WHERE sr.coord_system_id = cs.coord_system_id
    AND sr.seq_region_id = sra.seq_region_id
    AND sra.attrib_type_id = at.attrib_type_id
    AND at.code = 'toplevel'
  );

  # ref
  $sth = $dbh->{'ref'}->prepare($sql);
  $sth->execute;
  while (my ($cs) = $sth->fetchrow_array) {
    $common_cs{$cs}++;
  }
  $sth->finish;
  
  # alt
  $sth = $dbh->{'alt'}->prepare($sql);
  $sth->execute;
  while (my ($cs) = $sth->fetchrow_array) {
    $common_cs{$cs}++;
  }
  $sth->finish;

  # now determine the common ones
  foreach my $cs (sort keys %common_cs) {
    push @coord_systems, $cs if ($common_cs{$cs} > 1);
  }
}

unless (@coord_systems) {
  $support->error("No common toplevel coord_systems found.\n");
}

my $cs_string = join("', '", @coord_systems);
$cs_string = "'$cs_string'";
$support->log("Will use these coord_systems: $cs_string\n", 1);

# determine max(seq_region_id) and max(coord_system_id) in Ensembl
$sql = qq(SELECT MAX(seq_region_id) FROM seq_region);
$sth = $dbh->{'ref'}->prepare($sql);
$sth->execute;
my ($max_sri) = $sth->fetchrow_array;
my $sri_adjust = 1+$max_sri;

$sql = qq(SELECT MAX(coord_system_id) FROM coord_system);
$sth = $dbh->{'ref'}->prepare($sql);
$sth->execute;
my ($max_csi) = $sth->fetchrow_array;
my $csi_adjust = 10**(length($max_csi));

my $ref_db = $support->param('dbname');
my $alt_assembly = $support->param('altassembly');

# fetch and insert alternative seq_regions with adjusted seq_region_id and
# coord_system_id
$sql = qq(
    INSERT IGNORE INTO $ref_db.seq_region
    SELECT
        sr.seq_region_id+$sri_adjust,
        sr.name,
        sr.coord_system_id+$csi_adjust,
        sr.length
    FROM seq_region sr, coord_system cs, seq_region_attrib sra, attrib_type at
    WHERE sr.coord_system_id = cs.coord_system_id
    AND sr.seq_region_id = sra.seq_region_id
    AND sra.attrib_type_id = at.attrib_type_id
    AND at.code = 'toplevel'
    AND cs.name IN ($cs_string)
    AND cs.version = '$alt_assembly'
);
my $c = $dbh->{'alt'}->do($sql);
$support->log_stamped("Done loading $c seq_regions.\n\n");

## Add in any codon table attributes as HCs complain otherwise
$sql = qq(
    INSERT IGNORE INTO $ref_db.seq_region_attrib
    SELECT
        sr.seq_region_id+$sri_adjust,
        codon_at.attrib_type_id,
        codon_sra.value
    FROM seq_region sr, coord_system cs, seq_region_attrib sra, attrib_type at, seq_region_attrib codon_sra, attrib_type codon_at
    WHERE sr.coord_system_id = cs.coord_system_id
    AND sr.seq_region_id = sra.seq_region_id
    AND sra.attrib_type_id = at.attrib_type_id
    AND sr.seq_region_id = codon_sra.seq_region_id
    AND codon_sra.attrib_type_id = codon_at.attrib_type_id
    AND at.code = 'toplevel'
    AND cs.name IN ($cs_string)
    AND cs.version = '$alt_assembly'
    AND codon_at.code = 'codon_table';
);
$c = 0;
$c = $dbh->{'alt'}->do($sql);
$support->log_stamped("Done loading $c seq_region codon_table attributes.\n\n");

#####
# add appropriate entries to coord_system
#
$support->log_stamped("Adding coord_system entries...\n");
$c = 0;
foreach my $cs (@coord_systems) {
  $sql = qq(
      INSERT IGNORE INTO $ref_db.coord_system
      SELECT coord_system_id+$csi_adjust, 1,name, version,
        (SELECT MAX(rank)+1 FROM $ref_db.coord_system), ''
      FROM coord_system
      WHERE name = '$cs'
      AND version = '$alt_assembly'
  );
  $c += $dbh->{'alt'}->do($sql);
}
$support->log_stamped("Done adding $c coord_system entries.\n\n");

#####
# add assembly.mapping to meta table
#
$support->log_stamped("Adding assembly.mapping entry to meta table...\n");
$c = 0;
foreach my $cs (@coord_systems) {
  my $mappingstring = "$cs:".$support->param('assembly').
      "#$cs:".$support->param('altassembly');
  $sql = qq(
      INSERT IGNORE INTO meta (meta_key, meta_value)
      VALUES ('assembly.mapping', '$mappingstring')
  );
  $c += $dbh->{'ref'}->do($sql);
}
$support->log_stamped("Done inserting $c meta entries.\n\n");

# finish logfile
$support->finish_log;

