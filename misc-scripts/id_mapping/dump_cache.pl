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
# Change the PATH in the myRun.ksh script if you want to use another perl.

=head1 NAME


=head1 SYNOPSIS

.pl [arguments]

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

  -i, --interactive=0|1               run script interactively (default: true)
  -n, --dry_run, --dry=0|1            don't write results to database
  -h, --help, -?                      print help (this message)

=head1 DESCRIPTION

Use --sourceschema and --targetschema to specify a schema version (default:
latest). This will be used to determine the subroutine to build the cache. By
default, &build_cache_latest() is run which uses Bio::EnsEMBL::IdMapping::Cache
to read from the database and write the cache.  An alternative subroutine can
use a different module for that, which will usually inherit from the former and
overwrite Cache->build_cache(). This is useful for backwards compatibility with
older schema versions. Once the cache is built, no API access is needed,
therefore the ID mapping application is independent of the underlying database
schema.



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
use Bio::EnsEMBL::Utils::ScriptUtils qw(inject path_append);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

# parse configuration and commandline arguments
my $conf = new Bio::EnsEMBL::Utils::ConfParser(
  -SERVERROOT => "$Bin/../../..",
  -DEFAULT_CONF => "$Bin/default.conf"
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

# set default logpath
unless ($conf->param('logpath')) {
  $conf->param('logpath', path_append($conf->param('basedir'), 'log'));
}

# get log filehandle and print heading and parameters to logfile
my $logger = new Bio::EnsEMBL::Utils::Logger(
  -LOGFILE      => $conf->param('logfile'),
  -LOGAUTO      => $conf->param('logauto'),
  -LOGAUTOBASE  => 'dump_cache',
  -LOGAUTOID    => $conf->param('logautoid'),
  -LOGPATH      => $conf->param('logpath'),
  -LOGAPPEND    => $conf->param('logappend'),
  -LOGLEVEL     => $conf->param('loglevel'),
  -IS_COMPONENT => $conf->param('is_component'),
);

# initialise log
$logger->init_log($conf->list_param_values);

# determin cache method to use.
# this can be used to support different caching strategies or access to old
# database schemas.
my $cache_method = $conf->param('cache_method') || 'build_cache_auto';
no strict 'refs';
my $retval = &$cache_method;

# finish logfile
$logger->finish_log;

exit($retval);


### END main ###


sub build_cache_auto {
  # load the cache implementation
  my $cache_impl = 'Bio::EnsEMBL::IdMapping::Cache';
  inject($cache_impl);

  my $cache = $cache_impl->new(
    -LOGGER       => $logger,
    -CONF         => $conf,
  );

  $logger->debug("\nChecking number of toplevel seq_regions...\n");
  my $max = 0;

  foreach my $dbtype (qw(source target)) {
    my $num = scalar(@{ $cache->slice_names($dbtype) });
    $max = $num if ($num > $max);
    $logger->debug("$dbtype: $num.\n", 1);
  }

  my $threshold = $conf->param('build_cache_auto_threshold') || 100;
  my $retval;

  if ($max > $threshold) {
    $logger->debug("\nWill use build_cache_all.\n");
    $retval = &build_cache_all;
  } else {
    $logger->debug("\nWill use build_cache_by_seq_region.\n");
    $retval = &build_cache_by_seq_region;
  }

  return $retval;
}


sub build_cache_by_seq_region {

  my %jobs = ();

  # create empty directory for logs
  my $logpath = path_append($conf->param('logpath'), 'dump_by_seq_region');
  system("rm -rf $logpath") == 0 or
    $logger->error("Unable to delete lsf log dir $logpath: $!\n");
  system("mkdir -p $logpath") == 0 or
    $logger->error("Can't create lsf log dir $logpath: $!\n");

  # load the cache implementation
  my $cache_impl = 'Bio::EnsEMBL::IdMapping::Cache';
  inject($cache_impl);

  my $cache = $cache_impl->new(
    -LOGGER       => $logger,
    -CONF         => $conf,
  );

  # submit jobs to lsf
  foreach my $dbtype (qw(source target)) {

    $logger->info("\n".ucfirst($dbtype)." db...\n", 0, 'stamped');

    # determine which slices need to be done
    my $filename = "$dbtype.dump_cache.slices.txt";
    open(my $fh, '>', "$logpath/$filename") or
      throw("Unable to open $logpath/$filename for writing: $!");
    
    my $num_jobs = 0;

    foreach my $slice_name (@{ $cache->slice_names($dbtype) }) {
      my $type = "$dbtype.$slice_name";
      unless ($cache->cache_file_exists($type)) {
        print $fh "$slice_name\n";
        $num_jobs++;
      }
    }

    close($fh);

    unless ($num_jobs) {
      $logger->info("All cache files for $dbtype exist.\n");
      next;
    }

    # build lsf command
    my $lsf_name = 'dump_by_seq_region_'.time;
    my $concurrent = $conf->param('build_cache_concurrent_jobs') || 200;

    my $options = $conf->create_commandline_options(
        logauto       => 1,
        logautobase   => "dump_by_seq_region",
        interactive   => 0,
        is_component  => 1,
        dbtype        => $dbtype,
        cache_impl    => $cache_impl,
    );

    my $cmd = qq{./dump_by_seq_region.pl $options --index \$LSB_JOBINDEX};

    my $pipe =
        qq{|bsub -J '$lsf_name\[1-$num_jobs\]\%$concurrent' }
      . qq{-o $logpath/dump_by_seq_region.$dbtype.\%I.out }
      . qq{-e $logpath/dump_by_seq_region.$dbtype.\%I.err }
      . $conf->param('lsf_opt_dump_cache');

    # run lsf job array
    $logger->info("\nSubmitting $num_jobs jobs to lsf.\n");
    $logger->debug("$cmd\n\n");
    $logger->debug("$pipe\n\n");

    local *BSUB;
    open BSUB, $pipe or
      $logger->error("Could not open open pipe to bsub: $!\n");

    print BSUB $cmd;
    $logger->error("Error submitting jobs: $!\n")
      unless ($? == 0); 
    close BSUB;

    # submit dependent job to monitor finishing of jobs
    $logger->info("Waiting for jobs to finish...\n", 0, 'stamped');

    my $dependent_job =
      qq{bsub -K -w "ended($lsf_name)" -q production-rh7 } .
      qq{-M 100 -R 'select[mem>100]' -R 'rusage[mem=100]' } .
      qq{-o $logpath/dump_cache.$dbtype.depend.out /bin/true};

    system($dependent_job) == 0 or
      $logger->error("Error submitting dependent job: $!\n");

    $logger->info("All jobs finished.\n", 0, 'stamped');

    # check for lsf errors
    sleep(5);
    my $err;
    foreach my $i (1..$num_jobs) {
      $err++ unless (-e "$logpath/dump_by_seq_region.$dbtype.$i.success");
    }

    if ($err) {
      $logger->error("At least one of your jobs failed.\nPlease check the logfiles at $logpath for errors.\n");
      return 1;
    }

  }

  return 0;
}


sub build_cache_all {

  # load the cache implementation
  my $cache_impl = 'Bio::EnsEMBL::IdMapping::Cache';
  inject($cache_impl);

  my $cache = $cache_impl->new(
    -LOGGER       => $logger,
    -CONF         => $conf,
  );

  # submit jobs to lsf
  foreach my $dbtype (qw(source target)) {

    $logger->info("\n".ucfirst($dbtype)." db...\n", 0, 'stamped');
    $logger->info("Building cache for whole genome...\n");

    my $i = 0;
    my $size = 0;
    ($i, $size) = $cache->build_cache_all($dbtype);

    $logger->info("Done with $dbtype (genes: $i, filesize: $size).\n", 0,
      'stamped');
  }

  return 0;
}


