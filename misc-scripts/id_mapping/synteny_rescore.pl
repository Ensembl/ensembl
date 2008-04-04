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
use Bio::EnsEMBL::IdMapping::Cache;
use Bio::EnsEMBL::IdMapping::SyntenyFramework;
use Bio::EnsEMBL::IdMapping::ScoredMappingMatrix;

# parse configuration and commandline arguments
my $conf = new Bio::EnsEMBL::Utils::ConfParser(
  -SERVERROOT => "$Bin/../../..",
  -DEFAULT_CONF => "$Bin/default.conf"
);

$conf->parse_options(
  'dumppath|dump_path=s' => 1,
  'index|i=n' => 1,
);

# append job index to logfile name
my $index = $conf->param('index');
my $logautobase = ($conf->param('logautobase') || 'synteny_rescore_') . $index;

# get log filehandle and print heading and parameters to logfile
my $logger = new Bio::EnsEMBL::Utils::Logger(
  -LOGFILE      => $conf->param('logfile'),
  -LOGAUTO      => $conf->param('logauto'),
  -LOGAUTOBASE  => $logautobase,
  -LOGAUTOID    => $conf->param('logautoid'),
  -LOGPATH      => $conf->param('logpath'),
  -LOGAPPEND    => $conf->param('logappend'),
  -LOGLEVEL     => $conf->param('loglevel'),
  -IS_COMPONENT => 1,
);

# loading cache from file
my $cache = Bio::EnsEMBL::IdMapping::Cache->new(
  -LOGGER       => $logger,
  -CONF         => $conf,
);
$cache->read_instance_from_file;

# load SyntenyFramework and gene ScoredMappingMatrix from files
my $dump_path = $conf->param('dumppath');

my $sf = Bio::EnsEMBL::IdMapping::SyntenyFramework->new(
  -DUMP_PATH    => path_append($dump_path, 'mapping'),
  -CACHE_FILE   => 'synteny_framework.ser',
  -LOGGER       => $logger,
  -CONF         => $conf,
  -CACHE        => $cache,
);
$sf->read_from_file;

my $gene_scores = Bio::EnsEMBL::IdMapping::ScoredMappingMatrix->new(
  -DUMP_PATH   => path_append($dump_path, 'matrix/synteny_rescore'),
  -CACHE_FILE  => "gene_matrix_synteny$index.ser",
  -AUTO_LOAD   => 1,
);

# synteny rescore and serialise results to file
$gene_scores = $sf->rescore_gene_matrix($gene_scores);
$gene_scores->write_to_file;

# set flag to indicate everything went fine
my $success_file = $conf->param('logpath')."/synteny_rescore_$index.success";
open(TMPFILE, '>', $success_file) and close TMPFILE
  or die "Can't open $success_file for writing: $!";

# log success
$logger->finish_log;

