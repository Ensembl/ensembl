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
  'dumppath|dump_path=s' => 1,
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
);

# set default logpath
unless ($conf->param('logpath')) {
  $conf->param('logpath', path_append($conf->param('dumppath'), 'log'));
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

# if user wants to run via lsf, submit script with bsub (this will exit this
# instance of the script)
&bsubmit if ($conf->param('lsf'));

# initialise log
$logger->init_log($conf->list_param_values);

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

# run components, depending on mode
my $mode = $conf->param('mode') || 'normal';
my $sub = "run_$mode";
no strict 'refs';
&$sub;

# finish logfile
$logger->finish_log;

### END main ###


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
    logautoid => $logger->log_auto_id,
    interactive   => 0,
    lsf       => 0,
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

