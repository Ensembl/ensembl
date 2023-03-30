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

package Bio::EnsEMBL::IdMapping::Pipeline::InitCheck;

use strict;
use warnings;
use Data::Dumper;
use File::Basename;

use base qw/Bio::EnsEMBL::Hive::Process/;
use Bio::EnsEMBL::Utils::ConfParser;
use Bio::EnsEMBL::Utils::Logger;
use Bio::EnsEMBL::Utils::ScriptUtils qw(path_append);
use Bio::EnsEMBL::IdMapping::Cache;

my ($conf, $logger);

sub run {
  my ($self) = @_;

  my $config   = $self->param_required('config');
  my $mode     = $self->param_required('mode');
  my $logauto  = $self->param_required('logauto');
  my $no_check = $self->param_required('no_check');

  # Parse configuration and commandline arguments
  $conf = new Bio::EnsEMBL::Utils::ConfParser(
    -SERVERROOT => dirname($config)."/../../..",
    -DEFAULT_CONF => $config
  );

  $conf->parse_options(
    'sourcehost|source_host=s' => 1,
    'sourceport|source_port=n' => 1,
    'sourceuser|source_user=s' => 1,
    'sourcepass|source_pass=s' => 0,
    'sourcedbname|source_dbname=s' => 1,
    'targethost|target_host=s' => 1,
    'targetport|target_port=n' => 1,
    'targetuser|target_user=s' => 1,
    'targetpass|target_pass=s' => 0,
    'targetdbname|target_dbname=s' => 1,
    'mode=s' => 0,
    'basedir|basedir=s' => 1,
    'chromosomes|chr=s@' => 0,
    'region=s' => 0,
    'biotypes=s@' => 0,
    'biotypes_include=s@' => 0,
    'biotypes_exclude=s@' => 0,
    'cache_method=s' => 0,
    'build_cache_auto_threshold=n' => 0,
    'build_cache_concurrent_jobs=n' => 0,
    'min_exon_length|minexonlength=i' => 0,
    'exonerate_path|exoneratepath=s' => 1,
    'exonerate_threshold|exoneratethreshold=f' => 0,
    'exonerate_jobs|exoneratejobs=i' => 0,
    'exonerate_bytes_per_job|exoneratebytesperjob=f' => 0,
    'exonerate_extra_params|exonerateextraparams=s' => 0,
    'plugin_internal_id_mappers_gene=s@' => 0,
    'plugin_internal_id_mappers_transcript=s@' => 0,
    'plugin_internal_id_mappers_exon=s@' => 0,
    'mapping_types=s@' => 1,
    'plugin_stable_id_generator=s' => 0,
    'upload_events|uploadevents=s' => 0,
    'upload_stable_ids|uploadstableids=s' => 0,
    'upload_archive|uploadarchive=s' => 0,
    'lsf!' => 0,
    'lsf_opt_run|lsfoptrun=s' => 0,
    'lsf_opt_dump_cache|lsfoptdumpcache=s' => 0,
    'no_check!' => 0,
    'no_check_empty_tables' => 0,
  );

  $conf->param('mode', $mode);
  $conf->param('logauto', $logauto);
  $conf->param('no_check', $no_check);

  # Set default logpath
  unless ($conf->param('logpath')) {
    $conf->param('logpath', path_append($conf->param('basedir'), 'log'));
  }

  # Get log filehandle
  $logger = new Bio::EnsEMBL::Utils::Logger(
    -LOGFILE      => $conf->param('logfile'),
    -LOGAUTO      => $conf->param('logauto'),
    -LOGAUTOBASE  => 'run_all',
    -LOGPATH      => $conf->param('logpath'),
    -LOGAPPEND    => $conf->param('logappend'),
    -LOGLEVEL     => $conf->param('loglevel'),
  );

  # Initialise log
  $logger->init_log($conf->list_param_values);

  # Check configuration and resources
  unless ($conf->param('no_check')) {
    if ($self->init_check($mode) > 0) {
      $logger->error("Configuration check failed. See above for details.\n");
      $self->log()->error("Configuration check failed.");
    }

    if ($mode eq 'check_only') {
      $logger->info("Nothing else to do for 'check_only' mode. Exiting.\n");
      $self->log()->info("Nothing else to do for 'check_only' mode. Exiting.");
      $self->complete_early("Nothing else to do for 'check_only' mode. Exiting");
    }
  }

  my $logautoid = $logger->log_auto_id;
  $self->dataflow_output_id({
    config    => $config,
    logautoid => $logautoid,
    logauto   => $logauto,
    mode      => $mode
  }, 1);

  # Finish logfile
  $logger->finish_log;
}

sub init_check {
  my $self = shift;
  my $mode = shift;

  my $err = 0;
  my %valid_modes = (
    'check_only' => 1,
    'normal'     => 1,
    'upload'     => 1,
    'mapping'    => 1
  );

  $logger->info("Checking configuration...\n", 0, 'stamped');
  $self->log()->info("Checking configuration...");

  # Check for valid mode
  unless ($valid_modes{$mode}) {
    $logger->warning("Invalid mode: $mode.\n");
    $self->log()->warn("Invalid mode: $mode.");
    $err++;
  } else {
    $logger->debug("Run mode ok.\n");
    $self->log()->debug("Run mode ok.");
  }

  # Create the base directory, throw if this fails
  my $basedir = $conf->param('basedir');
  unless (-d $basedir) {
    if (system("mkdir -p $basedir") == 0) {
      $logger->debug("Base directory created successfully.\n");
      $self->log()->debug("Base directory created successfully.");
    } else {
      $logger->warning("Unable to create base directory $basedir: $!\n");
      $self->log()->warn("Unable to create base directory $basedir: $!");
      $err++;
    }
  }

  # Check db connection and permissions (SELECT for source, INSERT for target)
  my $cache = Bio::EnsEMBL::IdMapping::Cache->new(
    -LOGGER => $logger,
    -CONF   => $conf,
  );

  # Source db
  $err += $cache->check_db_connection('source');
  $err += $cache->check_db_read_permissions('source');
  
  # Target db
  $err += $cache->check_db_connection('target');
  $err += $cache->check_db_read_permissions('target');
  $err += $cache->check_db_write_permissions('target');
  
  # Check stable ID and archive tables in target db are empty
  $err += $cache->check_empty_tables('target');
  
  # Check both dbs have sequence
  $err += $cache->check_sequence('source');
  $err += $cache->check_sequence('target');
  
  # Check for required meta table entries
  $err += $cache->check_meta_entries('source');
  $err += $cache->check_meta_entries('target');

  $logger->info("Done.\n\n", 0, 'stamped');
  $self->log()->info("Done.");

  return $err;
}

1;
