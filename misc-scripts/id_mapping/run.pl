#!/software/bin/perl

=head1 NAME

run_all.pl - wrapper script to run the stable ID mapping

=head1 SYNOPSIS

run_all.pl [arguments]

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
use Bio::EnsEMBL::Utils::ScriptUtils qw(path_append);

my %valid_modes = (
  'check_only'  => 1,
  'normal'      => 1,
  'upload'      => 1,
);

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
  'mode=s' => 0,
  'basedir|basedir=s' => 1,
  'chromosomes|chr=s@' => 0,
  'region=s' => 0,
  'biotypes=s@' => 0,
  'min_exon_length|minexonlength=i' => 0,
  'exonerate_path|exoneratepath=s' => 1,
  'exonerate_threshold|exoneratethreshold=f' => 0,
  'exonerate_jobs|exoneratejobs=i' => 0,
  'exonerate_bytes_per_job|exoneratebytesperjob=f' => 0,
  'exonerate_extra_params|exonerateextraparams=s' => 0,
  'upload_events|uploadevents=s' => 0,
  'upload_stable_ids|uploadstableids=s' => 0,
  'upload_archive|uploadarchive=s' => 0,
  'lsf!' => 0,
  'lsf_opt_run|lsfoptrun=s' => 0,
  'lsf_opt_dump_cache|lsfoptdumpcache=s' => 0,
  'no_check!' => 0,
);

# set default logpath
unless ($conf->param('logpath')) {
  $conf->param('logpath', path_append($conf->param('basedir'), 'log'));
}

# get log filehandle and print heading and parameters to logfile
my $logger = new Bio::EnsEMBL::Utils::Logger(
  -LOGFILE      => $conf->param('logfile'),
  -LOGAUTO      => $conf->param('logauto'),
  -LOGAUTOBASE  => 'run_all',
  -LOGPATH      => $conf->param('logpath'),
  -LOGAPPEND    => $conf->param('logappend'),
  -LOGLEVEL     => $conf->param('loglevel'),
);

# initialise log
$logger->init_log($conf->list_param_values);

my $mode = $conf->param('mode') || 'normal';

# check configuration and resources.
# this is deliberately done before submitting to lsf (doesn't need much
# resources and you will know about config errors before waiting for job to
# run). the 'no_check' option prevents the checks to be re-run after automatic
# lsf submission
unless ($conf->param('no_check')) {
  unless (&init_check($mode)) {
    $logger->error("Configuration check failed. See above for details.\n");
  }

  if ($mode eq 'check_only') {
    $logger->info("Nothing else to do for 'check_only' mode. Exiting.\n");
    exit;
  }
}

# if user wants to run via lsf, submit script with bsub (this will exit this
# instance of the script)
&bsubmit if ($conf->param('lsf'));

# this script is only a wrapper and will run one or more components.
# define options for the components here.
my %options;
my $logautoid = $logger->log_auto_id;

$options{'dump_cache'} = $conf->create_commandline_options(
  logautoid     => $logautoid,
  logappend     => 1,
  interactive   => 0,
  is_component  => 1,
);

$options{'id_mapping'} = $conf->create_commandline_options(
  logautoid     => $logautoid,
  logappend     => 1,
  interactive   => 0,
  is_component  => 1,
);

# $logger->info("Nothing to do, just testing.\n");
# exit;

# run components, depending on mode
my $sub = "run_$mode";
no strict 'refs';
&$sub;

# finish logfile
$logger->finish_log;

### END main ###
  # add one more job to 


sub init_check {
  my $mode = shift;

  my $ok = 1;

  $logger->info("Checking configuration...\n", 0, 'stamped');

  # check for valid mode
  unless ($valid_modes{$mode}) {
    $logger->warning("Invalid mode: $mode.\n");
    $ok = 0;
  }

  # create the base directory, throw if this fails
  my $basedir = $conf->param('basedir');
  unless (-d $basedir) {
    system("mkdir -p $basedir") == 0 or
      die("Unable to create directory $basedir: $!\n");
  }
  
  $logger->info("Done.\n\n", 0, 'stamped');

  return $ok;
}


sub run_check_only {
  # this is a pseudo mode which does nothing; only init_check() will be executed
  return;
}


sub run_normal {
  
  # dump cache files (this is done for all modes)
  &run_component('dump_cache', $options{'dump_cache'}, 'building cache');

  # ID mapping
  &run_component('id_mapping', $options{'id_mapping'}, 'Id mapping');

  # QC
  #&run_component('qc', $options{'qc'}, 'QC');

}


sub run_upload {
  # upload table data files into db
  # (delegate to id_mapping.pl which will do the right thing based on --mode)
  &run_component('id_mapping', $options{'id_mapping'}, 'uploading tables');
}


sub run_component {
  my $basename = shift;
  my $options = shift;
  my $logtext = shift;

  my $cmd = "$basename.pl";
  $logtext ||= $cmd;

  $logger->info("----- $logtext -----\n", 0, 'stamped');
  
  if ($logger->logauto) {
    $logger->info("See ${basename}_".$logger->log_auto_id.".log for logs.\n", 1);
  } elsif ($logger->logfile) {
    $logger->info("See below for logs.\n", 1);
  }

  system("./$cmd $options") == 0
    or $logger->error("Error running $cmd. Please see the respective logfile for more information.\n");
  
  $logger->info("----- done with $logtext -----\n\n", 0, 'stamped');
}


sub bsubmit {
  #
  # build bsub commandline
  #

  # automatically create a filename for lsf output
  my $cmd = 'bsub -o '.$conf->param('logpath');
  $cmd .= '/lsf_'.$logger->log_auto_id.'.out';

  # add extra lsf options as configured by the user
  $cmd .= ' '.$conf->param('lsf_opt_run');

  # this script's name
  $cmd .= " $0";

  # options for this script
  my $options = $conf->create_commandline_options(
    logautoid   => $logger->log_auto_id,
    interactive => 0,
    lsf         => 0,
    no_check    => 1,
  );
  $cmd .= " $options";

  #
  # execute bsub
  #
  print "\nRe-executing via lsf:\n";
  print "$cmd\n\n";

  exec($cmd) or die "Could not exec $0 via lsf: $!\n";
  #exit;
}

