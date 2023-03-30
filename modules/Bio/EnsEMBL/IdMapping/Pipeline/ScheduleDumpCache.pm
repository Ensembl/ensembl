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

package Bio::EnsEMBL::IdMapping::Pipeline::ScheduleDumpCache;

use strict;
use warnings;
use Data::Dumper;
use File::Basename;

use base qw/Bio::EnsEMBL::Hive::Process/;
use Bio::EnsEMBL::Utils::ConfParser;
use Bio::EnsEMBL::Utils::Logger;
use Bio::EnsEMBL::IdMapping::Cache;
use Bio::EnsEMBL::Utils::ScriptUtils qw(path_append);

my ($conf, $logger);

sub run {
  my ($self) = @_;

  my $config       = $self->param_required('config');
  my $logautoid    = $self->param('logautoid');
  my $logauto      = $self->param_required('logauto');
  my $cache_method = $self->param_required('cache_method');

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

  $conf->param('cache_method', $cache_method);
  $conf->param('logauto', $logauto);

  # Set default logpath
  unless ($conf->param('logpath')) {
    $conf->param('logpath', path_append($conf->param('basedir'), 'log'));
  }

  $logger = new Bio::EnsEMBL::Utils::Logger(
    -LOGFILE      => $conf->param('logfile'),
    -LOGAUTO      => $conf->param('logauto'),
    -LOGAUTOBASE  => 'dump_cache',
    -LOGAUTOID    => $logautoid,
    -LOGPATH      => $conf->param('logpath'),
    -LOGAPPEND    => $conf->param('logappend'),
    -LOGLEVEL     => $conf->param('loglevel'),
    -IS_COMPONENT => 1,
  );

  # Initialise log
  $logger->init_log($conf->list_param_values);

  my $retval;
  if ($cache_method eq 'build_cache_auto') {
    $retval = $self->build_cache_auto();
  } elsif ($cache_method eq 'build_cache_by_seq_region') {
    $retval = $self->build_cache_by_seq_region();
  } elsif ($cache_method eq 'build_cache_all') {
    $retval = $self->build_cache_all();
  } else {
    $logger->error("Invalid cache_method $cache_method.\n");
    $self->log()->error("Invalid cache_method $cache_method.");
  }

  if ($retval && scalar(@{$retval}) > 0) {
    foreach my $params (@{$retval}) {
      $self->dataflow_output_id($params, 2);
    }
  }

  # Finish logfile
  $logger->finish_log;
}

sub build_cache_auto {
  my $self = shift;

  my $max = 0;
  my $retval;

  # Load the cache implementation
  my $cache = Bio::EnsEMBL::IdMapping::Cache->new(
    -LOGGER => $logger,
    -CONF   => $conf,
  );

  $logger->debug("\nChecking number of toplevel seq_regions...\n");
  $self->log()->debug("Checking number of toplevel seq_regions...");

  foreach my $dbtype (qw(source target)) {
    my $num = scalar(@{ $cache->slice_names($dbtype) });
    $max = $num if ($num > $max);
    $logger->debug("$dbtype: $num.\n", 1);
    $self->log()->debug("$dbtype: $num.");
  }

  my $threshold = $conf->param('build_cache_auto_threshold') || 100;

  if ($max > $threshold) {
    $logger->debug("\nWill use build_cache_all.\n");
    $self->log()->debug("Will use build_cache_all.");
    $retval = $self->build_cache_all();
  } else {
    $logger->debug("\nWill use build_cache_by_seq_region.\n");
    $self->log()->debug("Will use build_cache_by_seq_region.");
    $retval = $self->build_cache_by_seq_region();
  }

  return $retval;
}


sub build_cache_by_seq_region {
  my $self = shift;

  my %jobs = ();
  my $retval;

  # Create empty directory for logs
  my $logpath = path_append($conf->param('logpath'), 'dump_by_seq_region');
  system("rm -rf $logpath") == 0 or
    $logger->error("Unable to delete lsf log dir $logpath: $!\n");
  system("mkdir -p $logpath") == 0 or
    $logger->error("Can't create lsf log dir $logpath: $!\n");

  # Load the cache implementation
  my $cache = Bio::EnsEMBL::IdMapping::Cache->new(
    -LOGGER => $logger,
    -CONF   => $conf,
  );

  # Create jobs
  foreach my $dbtype (qw(source target)) {
    $logger->info("\n".ucfirst($dbtype)." db...\n", 0, 'stamped');
    $self->log()->info(ucfirst($dbtype)." db...");

    # Determine which slices need to be done
    my $num_jobs = 0;
    my $dataflow_params;

    foreach my $slice_name (@{ $cache->slice_names($dbtype) }) {
      my $type = "$dbtype.$slice_name";
      unless ($cache->cache_file_exists($type)) {
        $num_jobs++;

        $dataflow_params = {
          dbtype  => $dbtype,
          slice_name => $slice_name
        };
        push(@{$retval}, $dataflow_params);
      }
    }

    unless ($num_jobs) {
      $logger->info("All cache files for $dbtype exist.\n");
      $self->log()->info("All cache files for $dbtype exist.");
      next;
    }
  }

  return $retval;
}


sub build_cache_all {
  my $self = shift;

  # Load the cache implementation
  my $cache = Bio::EnsEMBL::IdMapping::Cache->new(
    -LOGGER => $logger,
    -CONF   => $conf,
  );

  foreach my $dbtype (qw(source target)) {
    $logger->info("\n".ucfirst($dbtype)." db...\n", 0, 'stamped');
    $logger->info("Building cache for whole genome...\n");
    $self->log()->info(ucfirst($dbtype)." db...");
    $self->log()->info("Building cache for whole genome...");

    my $i = 0;
    my $size = 0;
    ($i, $size) = $cache->build_cache_all($dbtype);

    $logger->info("Done with $dbtype (genes: $i, filesize: $size).\n", 0,
      'stamped');
    $self->log()->info("Done with $dbtype (genes: $i, filesize: $size).");
  }

  return 0;
}

1;
