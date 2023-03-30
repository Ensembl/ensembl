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

package Bio::EnsEMBL::IdMapping::Pipeline::DumpCacheBySeqRegion;

use strict;
use warnings;
use Data::Dumper;
use File::Basename;

use base qw/Bio::EnsEMBL::Hive::Process/;
use Bio::EnsEMBL::Utils::ConfParser;
use Bio::EnsEMBL::Utils::Logger;
use Bio::EnsEMBL::IdMapping::Cache;
use Bio::EnsEMBL::Utils::ScriptUtils qw(path_append);

sub run {
  my ($self) = @_;

  my $config     = $self->param_required('config');
  my $logautoid  = $self->param('logautoid');
  my $dbtype     = $self->param('dbtype');
  my $slice_name = $self->param('slice_name');
  my $logauto    = $self->param_required('logauto');

  # Parse configuration and commandline arguments
  my $conf = new Bio::EnsEMBL::Utils::ConfParser(
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
    'basedir|basedir=s' => 1,
    'chromosomes|chr=s@' => 0,
    'region=s' => 0,
    'biotypes=s@' => 0,
    'biotypes_include=s@' => 0,
    'biotypes_exclude=s@' => 0,
    'lsf_opt_dump_cache|lsfoptdumpcache=s' => 0,
    'cache_method=s' => 0,
    'build_cache_auto_threshold=n' => 0,
    'build_cache_concurrent_jobs=n' => 0,
  );

  $conf->param('logauto', $logauto);

  # Set default logpath
  unless ($conf->param('logpath')) {
    $conf->param('logpath', path_append($conf->param('basedir'), 'log'));
  }
  
  my $logger = new Bio::EnsEMBL::Utils::Logger(
    -LOGFILE      => $conf->param('logfile'),
    -LOGAUTO      => $conf->param('logauto'),
    -LOGAUTOBASE  => 'dump_by_seq_region',
    -LOGAUTOID    => $logautoid,
    -LOGPATH      => $conf->param('logpath'),
    -LOGAPPEND    => $conf->param('logappend'),
    -LOGLEVEL     => $conf->param('loglevel'),
    -IS_COMPONENT => 1,
  );

  # Initialise log
  $logger->init_log($conf->list_param_values);

  my $cache = Bio::EnsEMBL::IdMapping::Cache->new(
    -LOGGER => $logger,
    -CONF   => $conf,
  );

  # Build the cache for this slice
  $cache->build_cache_by_slice($dbtype, $slice_name);

  # Finish logfile
  $logger->finish_log;
}

1;
