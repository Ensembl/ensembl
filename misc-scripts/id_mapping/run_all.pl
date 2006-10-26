#!/usr/local/bin/perl

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

$| = 1;

my $conf = new Bio::EnsEMBL::Utils::ConfParser(
  -SERVERROOT => $SERVERROOT,
);

# parse options
$conf->param('default_conf', './default.conf');
$conf->parse_common_options(@_);
$conf->parse_extra_options(qw(
  sourcehost|source_host=s
  sourceport|source_port=n
  sourceuser|source_user=s
  sourcepass|source_pass=s
  sourcedbname|source_dbname=s
  targethost|target_host=s
  targetport|target_port=n
  targetuser|target_user=s
  targetpass|target_pass=s
  targetdbname|target_dbname=s
  mode=s
  dumppath|dump_path=s
  cachefile|cache_file=s
  chromosomes|chr=s@
  region=s
  biotypes=s@
  min_exon_length|minexonlength=i
  exonerate_path|exoneratepath=s
  exonerate_threshold|exoneratethreshold=i
  exonerate_jobs|exoneratejobs=i
  exonerate_bytes_per_job|exoneratebytesperjob=i
));
$conf->allowed_params(
  $conf->get_common_params,
  qw(
    sourcehost sourceport sourceuser sourcepass sourcedbname
    targethost targetport targetuser targetpass targetdbname
    mode
    dumppath cachefile
    chromosomes region biotypes
    min_exon_length
    exonerate_path exonerate_threshold exonerate_jobs exonerate_byte_per_job
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
  -LOGFILE    => $conf->param('logfile'),
  -LOGPATH    => $conf->param('logpath'),
  -LOGAPPEND  => $conf->param('logappend'),
  -VERBOSE    => $conf->param('verbose'),
);

# initialise log
$logger->init_log($conf->list_all_params);

# check required parameters were set
$conf->check_required_params(
  qw(
    sourcehost sourceport sourceuser sourcedbname
    targethost targetport targetuser targetdbname
    dumppath
  )
);

# define options for the components
my %options;

$options{'dump_cache'} = $conf->create_commandline_options(
    -ALLOWED_PARAMS => 1,
    -REPLACE => {
        interactive => 0,
        is_component => 1,
    },
    -EXCLUDE => [qw(mode min_exon_length exonerate_path exonerate_threshold
                    exonerate_jobs exonerate_byte_per_job)]
);

$options{'id_mapping'} = $conf->create_commandline_options(
    -ALLOWED_PARAMS => 1,
    -REPLACE => {
        interactive => 0,
        is_component => 1,
    },
    -EXCLUDE => [qw(
        sourcehost sourceport sourceuser sourcepass sourcedbname
        targethost targetport targetuser targetpass targetdbname
    )]
);

# dump cache files (this is done for all modes)
&run_component('dump_cache.pl', $options{'dump_cache'}, 'building cache');

# run components, depending on mode
my $mode = $conf->param('mode') || 'normal';
my $sub = "run_$mode";
no strict 'refs';
&$sub;


# finish logfile
$logger->finish_log;


### END main ###


sub run_normal {
  # ID mapping
  &run_component('id_mapping.pl', $options{'id_mapping'}, 'Id mapping');

  # archive
  #&run_component('archive.pl', $options{'archive'}, 'archive');

  # reporting
  #&run_component('report.pl', $options{'report'}, 'creating report');

  # QC
  #&run_component('qc.pl', $options{'qc'}, 'QC');
}


sub run_component {
  my $cmd = shift;
  my $options = shift;
  my $logtext = shift || $cmd;

  $logger->log_stamped("----- $logtext -----\n");
  
  system("./$cmd $options") == 0
    or $logger->log_error("Error running $cmd. Please see the respective logfile for more information.\n");
  
  $logger->log_stamped("----- done with $logtext -----\n\n");
}

