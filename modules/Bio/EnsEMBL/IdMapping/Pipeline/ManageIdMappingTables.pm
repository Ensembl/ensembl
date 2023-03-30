=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2023] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

package Bio::EnsEMBL::IdMapping::Pipeline::ManageIdMappingTables;

use strict;
use warnings;
no warnings 'uninitialized';
use Data::Dumper;
use File::Basename;

use base qw/Bio::EnsEMBL::Hive::Process/;
use Bio::EnsEMBL::Utils::ConfParser;
use Bio::EnsEMBL::Utils::Logger;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::ScriptUtils qw(path_append);

my %suffnum = ();
my ($sfx, $dbh);
my ($conf, $logger);

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

sub run {
  my ($self) = @_;

  my $config             = $self->param_required('config');
  my $drop_backup_tables = $self->param('drop_backup_tables');
  my $backup_tables      = $self->param('backup_tables');
  my $delete_from_tables = $self->param('delete_from_tables');

  # Parse configuration and commandline arguments
  $conf = new Bio::EnsEMBL::Utils::ConfParser(
    -SERVERROOT => dirname($config)."/../../..",
    -DEFAULT_CONF => $config
  );

  $conf->parse_options(
    'targethost|target_host=s' => 1,
    'targetport|target_port=n' => 1,
    'targetuser|target_user=s' => 1,
    'targetpass|target_pass=s' => 0,
    'targetdbname|target_dbname=s' => 1,
  );

  # Set default logpath
  unless ($conf->param('logpath')) {
    $conf->param('logpath', path_append($conf->param('basedir'), 'log'));
  }

  # Get log filehandle
  $logger = new Bio::EnsEMBL::Utils::Logger(
    -LOGFILE   => $conf->param('logfile'),
    -LOGPATH   => $conf->param('logpath'),
    -LOGAPPEND => $conf->param('logappend'),
    -VERBOSE   => $conf->param('verbose'),
  );

  # Initialise log
  $logger->init_log($conf->list_param_values);

  # Connect to database and get adaptors
  my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host   => $conf->param('targethost'),
    -port   => $conf->param('targetport'),
    -user   => $conf->param('targetuser'),
    -pass   => $conf->param('targetpass'),
    -dbname => $conf->param('targetdbname'),
    -group  => 'core',
  );
  $dba->dnadb($dba);

  $dbh = $dba->dbc->db_handle;

  # First check which tables are populated
  $self->list_table_counts();

  # Then look for existing backup tables
  $sfx = $self->list_backup_counts();

  # Drop backup tables
  if ($drop_backup_tables) {
    $self->drop_backup_tables();
  }

  # Backup current tables
  if ($backup_tables) {
    $self->backup_tables($sfx);
  }

  # Delete from tables
  if ($delete_from_tables) {
    $self->delete_from_tables();
  }

  # Finish logfile
  $logger->finish_log;
}

sub list_table_counts {
  my $self = shift;

  $logger->info("Current table counts:\n\n");
  $self->log()->info("Current table counts:");
  $self->list_counts( [@tables] );
}

sub list_backup_counts {
  my $self = shift;
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
  $self->log()->info("Backup tables found:") if (%suffnum);

  foreach my $num ( sort keys %suffnum ) {
    my @t = ();

    foreach my $table (@tables) {
      my $thetable = $table;
      if ( $table =~ /^([^_]+)_stable_id/ ) {
        $thetable = $1;
      }
      push @t, "${thetable}_bak_$num";
    }

    $self->list_counts( [@t]);
    $logger->info("\n");
  }

  my $sfx = '_bak_' . ++$new_num;
  return $sfx;
}

sub list_counts {
  my $self = shift;
  my $tabs = shift;

  unless ( $tabs and ref($tabs) eq 'ARRAY' ) {
    throw("Need an arrayref.");
  }

  $logger->info( sprintf( "%-30s%-8s\n", qw(TABLE COUNT) ), 1 );
  $logger->info( ( '-' x 38 ) . "\n", 1 );
  $self->log()->info(sprintf( "%-30s%-8s\n", qw(TABLE COUNT) ));
  $self->log()->info(( '-' x 38 ));

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
    $self->log()->info(sprintf( $fmt, $thetable, $count ));
  }

  $logger->info("\n");
}

sub drop_backup_tables {
  my $self = shift;

  foreach my $num ( sort keys %suffnum ) {
    my $suffix = "_bak_$num";

    foreach my $table (@tables) {
      my $thetable = $table;
      if ( $table =~ /^([^_]+)_stable_id/ ) {
        $thetable = $1;
      }
      my $bak_table = "${thetable}${suffix}";
      $logger->info( "$bak_table\n", 1 );
      $self->log()->info($bak_table);
      unless ( $conf->param('dry_run') ) {
        $dbh->do(qq(DROP TABLE $bak_table));
      }
    }

    # Remove the suffix number
    delete $suffnum{$num};
  }

  $logger->info("\n");

  # Recalculate the suffix number to use for current backup
  my $max_num = reverse sort keys %suffnum;
  $sfx = '_bak_' . ++$max_num;
}

sub backup_tables {
  my $self = shift;
  my $sfx = shift;

  throw("Need a backup table suffix.") unless ( defined($sfx) );

  $logger->info(qq(\nWill use '$sfx' as suffix for backup tables\n));
  $self->log()->info("Will use '$sfx' as suffix for backup tables");

  $logger->info(qq(\nBacking up tables...\n));
  $self->log()->info("Backing up tables...");

  my $fmt1 = "%-30s";
  my $fmt2 = "%8d\n";

  foreach my $table (@tables) {
    my $thetable = $table;
    if ( $table =~ /^([^_]+)_stable_id/ ) {
      $thetable = $1;
    }
    $logger->info( sprintf( $fmt1, $thetable ), 1 );
    $self->log()->info(sprintf( $fmt1, $thetable ));
    my $c = 0;
    if ( !$conf->param('dry_run') &&
         $dbh->do(qq(CREATE TABLE ${thetable}${sfx} LIKE ${thetable})) )
    {
      $c = $dbh->do(
           qq(INSERT INTO ${thetable}${sfx} SELECT * FROM ${thetable}));
    }
    $logger->info( sprintf( $fmt2, $c ) );
    $self->log()->info(sprintf( $fmt2, $c ));
  }

  $logger->info(qq(Done.\n));
  $self->log()->info("Done");
}

sub delete_from_tables {
  my $self = shift;
  my $fmt1 = "%-30s";
  my $fmt2 = "%8d\n";

  $logger->info(qq(\nDeleting from current tables...\n));
  $self->log()->info("Deleting from current tables...");

  foreach my $table (@tables) {
    my $thetable = $table;
    if ( $table =~ /^([^_]+)_stable_id/ ) {
      $thetable = $1;
    }
    $logger->info( sprintf( $fmt1, $thetable ), 1 );
    $self->log()->info(sprintf( $fmt1, $thetable ));
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
    $self->log()->info(sprintf( $fmt2, $c ));
  }

  $logger->info(qq(Done.\n));
  $self->log()->info("Done");
}

1;
