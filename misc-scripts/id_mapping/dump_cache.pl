#!/usr/local/ensembl/bin/perl

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

  -v, --verbose=0|1                   verbose logging (default: false)
  -i, --interactive=0|1               run script interactively (default: true)
  -n, --dry_run, --dry=0|1            don't write results to database
  -h, --help, -?                      print help (this message)

=head1 DESCRIPTION

Use --oldschema and --newschema to specify a schema version (default: latest).
This will be used to determine the subroutine to build the cache. By default,
&build_cache_latest() is run which uses Bio::EnsEMBL::IdMapping::Cache to read
from the database and write the cache.  An alternative subroutine can use a
different module for that, which will usually inherit from the former and
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
use vars qw($SERVERROOT);

BEGIN {
    $SERVERROOT = "$Bin/../../..";
    unshift(@INC, "$SERVERROOT/ensembl/modules");
    unshift(@INC, "$SERVERROOT/bioperl-live");
}

use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Utils::ConfParser;
use Bio::EnsEMBL::Utils::Logger;
use Bio::EnsEMBL::Utils::ScriptUtils qw(dynamic_use);

$| = 1;

my $conf = new Bio::EnsEMBL::Utils::ConfParser(
  -SERVERROOT => $SERVERROOT,
);

# parse options
$conf->param('default_conf', './default.conf');
$conf->parse_common_options(@_);
$conf->parse_extra_options(qw(
  oldhost|old_host=s
  oldport|old_port=n
  olduser|old_user=s
  oldpass|old_pass=s
  olddbname|old_dbname=s
  newhost|new_host=s
  newport|new_port=n
  newuser|new_user=s
  newpass|new_pass=s
  newdbname|new_dbname=s
  dumppath|dump_path=s
  chromosomes|chr=s@
  region=s
  biotypes=s@
));
$conf->allowed_params(
  $conf->get_common_params,
  qw(
    oldhost oldport olduser oldpass olddbname
    newhost newport newuser newpass newdbname
    dumppath
    chromosomes region biotypes
  )
);

if ($conf->param('help') or $conf->error) {
    warn $conf->error if $conf->error;
    pod2usage(1);
}

# ask user to confirm parameters to proceed
$conf->confirm_params;

# get log filehandle and print heading and parameters to logfile
my $logger = new Bio::EnsEMBL::Utils::Logger(
  -LOGFILE      => $conf->param('logfile'),
  -LOGPATH      => $conf->param('logpath'),
  -LOGAPPEND    => $conf->param('logappend'),
  -VERBOSE      => $conf->param('verbose'),
  -IS_COMPONENT => $conf->param('is_component'),
);

# initialise log
$logger->init_log($conf->list_all_params);

# check required parameters were set
$conf->check_required_params(
  qw(
    oldhost oldport olduser olddbname
    newhost newport newuser newdbname
    dumppath
  )
);

my %jobs;

# create empty directory for logs
my $logpath = ($conf->param('logpath')||$conf->param('dumppath')).'/lsf';
system("rm -rf $logpath") == 0 or
  $logger->log_error("Unable to delete lsf log dir $logpath: $!\n");
system("mkdir $logpath") == 0 or
  $logger->log_error("Can't create lsf log dir $logpath: $!\n");

# submit jobs to lsf
foreach my $dbtype (qw(old new)) {
  
  $logger->log_stamped("\n".ucfirst($dbtype)." db...\n");

  my $schema = $conf->param("${dbtype}schema") || 'latest';
  my $cache_builder = "build_cache_$schema";
  no strict 'refs';
  &$cache_builder($dbtype);
}

# monitor progress
my $err;
my $total = scalar(keys %jobs);

while (keys %jobs) {
  foreach my $type (keys %jobs) {
    my $err_log = "$logpath/dump_by_seq_region.$type.err";

    # there's an error if the lsf error logfile has non-zero length
    $err++ if (-s $err_log);

    # the job has finished if you find the error logfile
    delete($jobs{$type}) if (-e $err_log);
  }

  $logger->log("Jobs waiting: ".scalar(keys %jobs)."/$total.\r");

  sleep(3) if (scalar(keys %jobs));
}

$logger->log("\n\n");

# check if anything went wrong
my $retval = 0;
if ($err) {
  $logger->log("At least one of your jobs failed.\n");
  $logger->log("Please check $logpath for errors.\n");
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
  
  $logger->log("Submitting jobs to lsf...\n", 1);

  foreach my $slice_name (@{ $cache->slice_names($dbtype) }) {
    &bsubmit($dbtype, $slice_name, $cache, $logpath);
  }

}


sub bsubmit {
  my $dbtype = shift;
  my $slice_name = shift;
  my $cache = shift;
  my $logpath = shift;
  
  $logger->log("$slice_name\n", 2);
  
  my $type = "$dbtype.$slice_name";

  unless ($cache->all_cache_files_exist($type)) {
    my $options = $conf->create_commandline_options(
        -ALLOWED_PARAMS => 1,
        -REPLACE => {
            interactive   => 0,
            is_component  => 1,
            dbtype        => $dbtype,
            slice_name    => $slice_name,
            cache_impl    => ref($cache),
        },
        -EXCLUDE => [qw(region chromosomes)]
    );
    
    my $cmd = "bsub ".
                "-o $logpath/dump_by_seq_region.$type.out ".
                "-e $logpath/dump_by_seq_region.$type.err ".
                "./dump_by_seq_region.pl $options";

    system("$cmd") == 0
      or $logger->log_error("Error running dump_by_seq_region.pl: $!");

    # mark job as submitted
    $jobs{$type} = 1;
  } 
  
}


