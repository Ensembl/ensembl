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

package Bio::EnsEMBL::IdMapping::Pipeline::AnalyseResults;

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
use Bio::EnsEMBL::IdMapping::ResultAnalyser;

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
    -LOGAUTOBASE  => 'analyse_results',
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

  # Final stats and mapping summary
  $self->analyse_results();

  # Finish logfile
  $logger->finish_log;
}

sub analyse_results {
  my $self = shift;

  $logger->info("Analysing results...\n", 0, 'stamped');
  $self->log()->info("Analysing results...");
  
  # Get a result analyser
  my $analyser = Bio::EnsEMBL::IdMapping::ResultAnalyser->new(
    -LOGGER       => $logger,
    -CONF         => $conf,
    -CACHE        => $cache
  );

  # Analyse results
  my $gene_mappings;
  $analyser->analyse($gene_mappings, $stable_id_mapper->get_all_stable_id_events('similarity'));

  # Write results to file
  $analyser->write_results_to_file;

  # Create click lists
  $analyser->create_clicklist;
  
  $logger->info("Done.\n\n", 0, 'stamped');
  $self->log()->info("Done.");

  # Mapping summary
  $logger->info("Creating mapping summary...\n", 0, 'stamped');
  $self->log()->info("Creating mapping summary...");
  $analyser->create_mapping_summary;
  $logger->info("Done.\n", 0, 'stamped');
  $self->log()->info("Done.");
}

1;
