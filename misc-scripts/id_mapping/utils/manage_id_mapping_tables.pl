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

# Don't change the above line.
# Change the PATH in the myManager.ksh script if you want to use another perl.

=head1 NAME

manage_id_mapping_tables.pl - script to delete (and optionally backup) ID
mapping results

=head1 SYNOPSIS

manage_id_mapping_tables.pl [arguments]

Required arguments:

  --dbname, db_name=NAME              database name NAME
  --host, --dbhost, --db_host=HOST    database host HOST
  --port, --dbport, --db_port=PORT    database port PORT
  --user, --dbuser, --db_user=USER    database username USER
  --pass, --dbpass, --db_pass=PASS    database passwort PASS

Optional arguments:

  --conffile, --conf=FILE             read parameters from FILE
                                      (default: conf/Conversion.ini)

  --logfile, --log=FILE               log to FILE (default: *STDOUT)
  --logpath=PATH                      write logfile to PATH (default: .)
  --logappend, --log_append           append to logfile (default: truncate)
  --loglevel=LEVEL                    define log level (default: INFO)

  -i, --interactive                   run script interactively (default: true)
  -n, --dry_run, --dry                don't write results to database
  -h, --help, -?                      print help (this message)

=head1 DESCRIPTION

This script will delete stable ID mapping data from a database. The script is
intended to be run interactively (your configuration will be overridden).

The tables that will be emptied are:

  gene_stable_id
  transcript_stable_id
  translation_stable_id
  exon_stable_id
  mapping_session
  stable_id_event
  gene_archive
  peptide_archive

Optionally (by interactive selection), the current tables can be backed up.
Backkup tables will get suffices of _bak_0, _bak_1, etc. (where the correct
number is determined automatically from existing backup tables). There is also
an option to drop existing backup tables.

Deleting from the current tables can also be skipped, so effectively this
script can do three different things (or any combination of them), depending on
your answers in the interactive process:

  - drop existing backup tables
  - backup current tables
  - delete from current tables


=head1 AUTHOR

Patrick Meidl <meidl@ebi.ac.uk>, Ensembl core API team

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<http://lists.ensembl.org/mailman/listinfo/dev>

=cut

use strict;
use warnings;
no warnings 'uninitialized';

use FindBin qw($Bin);
use Bio::EnsEMBL::Utils::ConfParser;
use Bio::EnsEMBL::Utils::Logger;
use Bio::EnsEMBL::Utils::ScriptUtils qw(user_proceed);
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my @tables = qw(
  gene_stable_id
  transcript_stable_id
  translation_stable_id
  exon_stable_id
  mapping_session
  stable_id_event
  gene_archive
  peptide_archive
);

my %suffnum = ();

# parse configuration and commandline arguments
my $conf =
  new Bio::EnsEMBL::Utils::ConfParser(-SERVERROOT => "$Bin/../../../..",
                                      -DEFAULT_CONF => "" );

$conf->parse_options( 'host=s'   => 1,
                      'port=n'   => 1,
                      'user=s'   => 1,
                      'pass=s'   => 0,
                      'dbname=s' => 1, );

# get log filehandle and print heading and parameters to logfile
my $logger =
  new Bio::EnsEMBL::Utils::Logger(
                                -LOGFILE   => $conf->param('logfile'),
                                -LOGPATH   => $conf->param('logpath'),
                                -LOGAPPEND => $conf->param('logappend'),
                                -VERBOSE   => $conf->param('verbose'),
  );

# always run interactively
$conf->param( 'interactive', 1 );

# initialise log
$logger->init_log( $conf->list_param_values );

# connect to database and get adaptors
my $dba =
  new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $conf->param('host'),
                                      -port   => $conf->param('port'),
                                      -user   => $conf->param('user'),
                                      -pass   => $conf->param('pass'),
                                      -dbname => $conf->param('dbname'),
                                      -group  => 'core', );
$dba->dnadb($dba);

my $dbh = $dba->dbc->db_handle;

# first check which tables are populated
&list_table_counts;

# then look for existing backup tables
my $sfx = &list_backup_counts;

# aks user if he wants to drop backup tables
if ( %suffnum
     and
     user_proceed(
          "Drop any backup tables? (you will be able chose which ones)",
          $conf->param('interactive'), 'n' ) )
{
  &drop_backup_tables;
}

# ask user if current tables should be backed up
if (user_proceed( "Backup current tables?", $conf->param('interactive'),
                  'y' ) )
{
  &backup_tables($sfx);
}

# delete from tables
if ( user_proceed( "Delete from current tables?",
                   $conf->param('interactive'),
                   'n' ) )
{
  &delete_from_tables;
}

# finish logfile
$logger->finish_log;

### END main ###

sub list_table_counts {
  $logger->info("Current table counts:\n\n");
  &list_counts( [@tables] );
}

sub list_backup_counts {
  my $new_num = -1;

  foreach my $table (@tables) {
    my $thetable = $table;
    if ( $table =~ /^([^_]+)_stable_id/ ) {
      $thetable = $1;
    }
    my $sth = $dbh->prepare(qq(SHOW TABLES LIKE "${thetable}_bak_%"));
    $sth->execute;

    while ( my ($bak) = $sth->fetchrow_array ) {
      if ($bak =~ /_bak_(\d+)$/) {
        my $num = $1;
        $suffnum{$num} = 1;

        $new_num = $num if ( $num > $new_num );
      }
    }

    $sth->finish;
  }

  $logger->info("Backup tables found:\n\n") if (%suffnum);

  foreach my $num ( sort keys %suffnum ) {
    my @t = ();

    foreach my $table (@tables) {
      my $thetable = $table;
      if ( $table =~ /^([^_]+)_stable_id/ ) {
        $thetable = $1;
      }
      push @t, "${thetable}_bak_$num";
    }

    &list_counts( [@t] );
    $logger->info("\n");
  }

  my $sfx = '_bak_' . ++$new_num;
  return $sfx;
} ## end sub list_backup_counts

sub list_counts {
  my $tabs = shift;

  unless ( $tabs and ref($tabs) eq 'ARRAY' ) {
    throw("Need an arrayref.");
  }

  $logger->info( sprintf( "%-30s%-8s\n", qw(TABLE COUNT) ), 1 );
  $logger->info( ( '-' x 38 ) . "\n", 1 );

  my $fmt = "%-30s%8d\n";

  foreach my $table (@$tabs) {
    my $sth;
    my $thetable = $table;
    if ( $table =~ /^([^_]+)_stable_id/ ) {
      $thetable = $1;
      $sth = $dbh->prepare(
        qq(SELECT COUNT(*) FROM $thetable WHERE stable_id IS NOT NULL));
    }
    else {
      $sth = $dbh->prepare(qq(SELECT COUNT(*) FROM $thetable));
    }
    $sth->execute;
    my $count = $sth->fetchrow_arrayref->[0];
    $sth->finish;

    $logger->info( sprintf( $fmt, $thetable, $count ), 1 );
  }

  $logger->info("\n");
} ## end sub list_counts

sub drop_backup_tables {

  foreach my $num ( sort keys %suffnum ) {
    my $suffix = "_bak_$num";
    if ( user_proceed( qq(Drop backup tables with suffix ${suffix}?),
                       $conf->param('interactive'), 'n' ) )
    {
      foreach my $table (@tables) {
        my $thetable = $table;
        if ( $table =~ /^([^_]+)_stable_id/ ) {
          $thetable = $1;
        }
        my $bak_table = "${thetable}${suffix}";
        $logger->info( "$bak_table\n", 1 );
        unless ( $conf->param('dry_run') ) {
          $dbh->do(qq(DROP TABLE $bak_table));
        }
      }

      # remove the suffix number
      delete $suffnum{$num};
    }
  }

  $logger->info("\n");

  # recalculate the suffix number to use for current backup
  my $max_num = reverse sort keys %suffnum;
  $sfx = '_bak_' . ++$max_num;
} ## end sub drop_backup_tables

sub backup_tables {
  my $sfx = shift;

  throw("Need a backup table suffix.") unless ( defined($sfx) );

  $logger->info(qq(\nWill use '$sfx' as suffix for backup tables\n));

  $logger->info(qq(\nBacking up tables...\n));

  my $fmt1 = "%-30s";
  my $fmt2 = "%8d\n";

  foreach my $table (@tables) {
    my $thetable = $table;
    if ( $table =~ /^([^_]+)_stable_id/ ) {
      $thetable = $1;
    }
    $logger->info( sprintf( $fmt1, $thetable ), 1 );
    my $c = 0;
    if ( !$conf->param('dry_run') &&
         $dbh->do(qq(CREATE TABLE ${thetable}${sfx} LIKE ${thetable})) )
    {
      $c = $dbh->do(
           qq(INSERT INTO ${thetable}${sfx} SELECT * FROM ${thetable}));
    }
    $logger->info( sprintf( $fmt2, $c ) );
  }

  $logger->info(qq(Done.\n));
} ## end sub backup_tables

sub delete_from_tables {
  my $fmt1 = "%-30s";
  my $fmt2 = "%8d\n";

  $logger->info(qq(\nDeleting from current tables...\n));

  foreach my $table (@tables) {
    my $thetable = $table;
    if ( $table =~ /^([^_]+)_stable_id/ ) {
      $thetable = $1;
    }
    $logger->info( sprintf( $fmt1, $thetable ), 1 );
    my $c = 0;
    unless ( $conf->param('dry_run') ) {
      if ( $table =~ /^([^_]+)_stable_id/ ) {
        $c = $dbh->do(qq(UPDATE $thetable SET stable_id=NULL));
      }
      else {
        $c = $dbh->do(qq(TRUNCATE $thetable));
      }
    }
    $logger->info( sprintf( $fmt2, $c ) );
  }

  $logger->info(qq(Done.\n));
} ## end sub delete_from_tables
