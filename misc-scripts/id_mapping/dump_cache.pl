#!/software/bin/perl

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


=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http://www.ensembl.org/info/about/code_licence.html for details.

=head1 AUTHOR

Patrick Meidl <meidl@ebi.ac.uk>, Ensembl core API team

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<ensembl-dev@ebi.ac.uk>

=cut

use strict;
use warnings;
no warnings 'uninitialized';

use FindBin qw($Bin);
use Bio::EnsEMBL::Utils::ConfParser;
use Bio::EnsEMBL::Utils::Logger;
use Bio::EnsEMBL::Utils::ScriptUtils qw(dynamic_use path_append);

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
  'lsf_opt_dump_cache|lsfoptdumpcache=s' => 0,
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

my %jobs;

# create empty directory for logs
my $logpath = path_append($conf->param('logpath'), 'dump_by_seq_region');
system("rm -rf $logpath") == 0 or
  $logger->error("Unable to delete lsf log dir $logpath: $!\n");
system("mkdir -p $logpath") == 0 or
  $logger->error("Can't create lsf log dir $logpath: $!\n");

# submit jobs to lsf
foreach my $dbtype (qw(source target)) {
  
  $logger->info("\n".ucfirst($dbtype)." db...\n", 0, 'stamped');

  my $schema = $conf->param("${dbtype}schema") || 'latest';
  my $cache_builder = "build_cache_$schema";
  no strict 'refs';
  &$cache_builder($dbtype);
}

# monitor progress
my $err;
my $total = scalar(keys %jobs);
my @types;

while (keys %jobs) {
  foreach my $type (keys %jobs) {
    my $err_log = "$logpath/dump_by_seq_region.$type.err";

    # there's an error if the lsf error logfile has non-zero length
    $err++ if (-s $err_log);

    # the job has finished if you find the error logfile
    delete($jobs{$type}) if (-e $err_log);
    push @types, $type;
  }

  $logger->info("Jobs waiting: ".scalar(keys %jobs)."/$total.\r");

  sleep(3) if (scalar(keys %jobs));
}

$logger->info("\n\n");

# check if anything went wrong
sleep(5);
foreach my $type (@types) {
  #warn "$logpath/dump_by_seq_region.$type.success\n";
  $err++ unless (-e "$logpath/dump_by_seq_region.$type.success");
}

my $retval = 0;
if ($err) {
  $logger->warning("At least one of your jobs failed.\n");
  $logger->warning("Please check $logpath and ".$logger->logpath."/".$logger->logfile." for errors.\n");
  $retval = 1;
}

# finish logfile
$logger->finish_log;

exit($retval);


### END main ###


sub build_cache_latest {
  my $dbtype = shift;

  my $cache_impl = 'Bio::EnsEMBL::IdMapping::Cache';

  dynamic_use($cache_impl);
  
  my $cache = $cache_impl->new(
    -LOGGER       => $logger,
    -CONF         => $conf,
  );

  my $dba = $cache->get_DBAdaptor($dbtype);
  my $sa = $dba->get_SliceAdaptor;
  
  $logger->info("Submitting jobs to lsf...\n", 1);

  foreach my $slice_name (@{ $cache->slice_names($dbtype) }) {
    &bsubmit($dbtype, $slice_name, $cache, $logpath);
  }

}


sub bsubmit {
  my $dbtype = shift;
  my $slice_name = shift;
  my $cache = shift;
  my $logpath = shift;
  
  $logger->info("$slice_name\n", 2);
  
  my $type = "$dbtype.$slice_name";

  unless ($cache->cache_file_exists($type)) {
    my $options = $conf->create_commandline_options(
        logauto       => 1,
        logautobase   => "dump_by_seq_region_${type}",
        interactive   => 0,
        is_component  => 1,
        dbtype        => $dbtype,
        slice_name    => $slice_name,
        cache_impl    => ref($cache),
    );

    my $cmd = "bsub ".
                "-o $logpath/dump_by_seq_region.$type.out ".
                "-e $logpath/dump_by_seq_region.$type.err ".
                $conf->param('lsf_opt_dump_cache') . " " .
                "./dump_by_seq_region.pl $options";

    system("$cmd") == 0
      or $logger->error("Error running dump_by_seq_region.pl: $!");

    # mark job as submitted
    $jobs{$type} = 1;
  } 
  
}


