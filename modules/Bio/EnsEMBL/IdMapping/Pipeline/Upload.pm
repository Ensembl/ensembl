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

package Bio::EnsEMBL::IdMapping::Pipeline::Upload;

use strict;
use warnings;
use Data::Dumper;
use File::Basename;

use base qw/Bio::EnsEMBL::Hive::Process/;
use Bio::EnsEMBL::Utils::ConfParser;
use Bio::EnsEMBL::Utils::Logger;
use Bio::EnsEMBL::IdMapping::Cache;
use Bio::EnsEMBL::Utils::ScriptUtils qw(path_append);
use Bio::EnsEMBL::IdMapping::StableIdMapper;

my ($conf, $logger);
my ($cache, $stable_id_mapper);

sub run {
  my ($self) = @_;

  my $config             = $self->param_required('config');
  my $logautoid          = $self->param('logautoid');
  my $logauto            = $self->param_required('logauto');

  # Parse configuration and commandline arguments
  $conf = new Bio::EnsEMBL::Utils::ConfParser(
    -SERVERROOT => dirname($config)."/../../..",
    -DEFAULT_CONF => $config
  );

  $conf->parse_options(
    'mode=s' => 0,
    'basedir|basedir=s' => 1,
    'chromosomes|chr=s@' => 0,
    'region=s' => 0,
    'biotypes=s@' => 0,
    'biotypes_include=s@' => 0,
    'biotypes_exclude=s@' => 0,
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
  );

  $conf->param('logauto', $logauto);

  # Set default logpath
  unless ($conf->param('logpath')) {
    $conf->param('logpath', path_append($conf->param('basedir'), 'log'));
  }
  
  $logger = new Bio::EnsEMBL::Utils::Logger(
    -LOGFILE      => $conf->param('logfile'),
    -LOGAUTO      => $conf->param('logauto'),
    -LOGAUTOBASE  => 'archive',
    -LOGAUTOID    => $logautoid,
    -LOGPATH      => $conf->param('logpath'),
    -LOGAPPEND    => $conf->param('logappend'),
    -LOGLEVEL     => $conf->param('loglevel'),
    -IS_COMPONENT => 1,
  );

  # Initialise log
  $logger->init_log($conf->list_param_values);

  # Loading cache from file
  $cache = Bio::EnsEMBL::IdMapping::Cache->new(
    -LOGGER         => $logger,
    -CONF           => $conf,
    -LOAD_INSTANCE  => 1,
  );

  # Get a stable ID mapper
  $stable_id_mapper = Bio::EnsEMBL::IdMapping::StableIdMapper->new(
    -LOGGER       => $logger,
    -CONF         => $conf,
    -CACHE        => $cache
  );

  # Upload table data files into db
  $self->upload_mapping_session_and_events();
  $self->upload_stable_ids();
  $self->upload_archive();

  # Finish logfile
  $logger->finish_log;
}

sub upload_mapping_session_and_events {
  my $self = shift;

  if ($conf->is_true('upload_events') && ! $conf->param('dry_run')) {
    $logger->info("Uploading mapping_session and stable_id_event tables...\n");
    $self->log()->info("Uploading mapping_session and stable_id_event tables...");

    my $i = 0;
    my $j = 0;
    
    $logger->info("mapping_session...\n", 1);
    $self->log()->info("mapping_session...");
    $i += $stable_id_mapper->upload_file_into_table('target', 'mapping_session', 'mapping_session.txt');
    $logger->info("$i\n", 1);
    $self->log()->info($i);
    
    $logger->info("stable_id_event...\n", 1);
    $self->log()->info("stable_id_event...");
    $j += $stable_id_mapper->upload_file_into_table('target', 'stable_id_event', 'stable_id_event_existing.txt');
    $j += $stable_id_mapper->upload_file_into_table('target', 'stable_id_event', 'stable_id_event_new.txt', 1);
    $j += $stable_id_mapper->upload_file_into_table('target', 'stable_id_event', 'stable_id_event_similarity.txt', 1);
    $logger->info("$j\n", 1);
    $self->log()->info($j);
    
    $logger->info("Done.\n\n");
    $self->log()->info("Done.");
  } else {
    $logger->info("Stable ID event and mapping session tables not uploaded.\n\n");
    $self->log()->info("Stable ID event and mapping session tables not uploaded.");
  }
}

sub upload_stable_ids {
  my $self = shift;

  if ($conf->is_true('upload_stable_ids') && ! $conf->param('dry_run')) {
    $logger->info("Uploading stable ID tables...\n");
    $self->log()->info("Uploading stable ID tables...");
    
    foreach my $t ($conf->param('mapping_types')) {
      $logger->info("${t}_stable_id...\n", 1);
      $self->log()->info("${t}_stable_id...");
      my $i = $stable_id_mapper->upload_file_into_table('target', "${t}_stable_id", "${t}_stable_id.txt");
      $logger->info("$i\n", 1);
      $self->log()->info($i);
    }
    
    $logger->info("Done.\n\n");
    $self->log()->info("Done.");
  } else {
    $logger->info("Stable ID tables not uploaded.\n\n");
    $self->log()->info("Stable ID tables not uploaded.");
  }
}

sub upload_archive {
  my $self = shift;

  if ($conf->is_true('upload_archive') && ! $conf->param('dry_run')) {
    $logger->info("Uploading gene and peptide tables...\n");
    $self->log()->info("Uploading gene and peptide tables...");
    
    foreach my $t (qw(gene peptide)) {
      $logger->info("${t}_archive...\n", 1);
      $self->log()->info("${t}_archive...");
      my $i = 0;
      $i += $stable_id_mapper->upload_file_into_table('target', "${t}_archive", "${t}_archive_existing.txt", 1);
      $i += $stable_id_mapper->upload_file_into_table('target', "${t}_archive", "${t}_archive_new.txt", 1);
      $logger->info("$i\n", 1);
      $self->log()->info($i);
    }
    
    $logger->info("Done.\n\n");
    $self->log()->info("Done.");
  } else {
    $logger->info("Gene and peptide archive tables not uploaded.\n\n");
    $self->log()->info("Gene and peptide archive tables not uploaded.");
  }
}

1;
