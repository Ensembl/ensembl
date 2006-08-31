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
use Bio::EnsEMBL::IdMapping::Cache;
use Bio::EnsEMBL::IdMapping::ExonScoreBuilder;

#use Devel::Size qw(size total_size);
#use Data::Dumper;
#$Data::Dumper::Indent = 1;


$| = 1;

my $conf = new Bio::EnsEMBL::Utils::ConfParser(
  -SERVERROOT => $SERVERROOT,
);

# parse options
$conf->param('default_conf', './default.conf');
$conf->parse_common_options(@_);
$conf->parse_extra_options(qw(
  mode=s
  dumppath|dump_path=s
  chromosomes|chr=s@
  region=s
  biotypes=s@
));
$conf->allowed_params(
  $conf->get_common_params,
  qw(
    mode
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
    dumppath
  )
);


# instance variables
my $exon_scores;
my $transcript_scores;
my $gene_scores;


# loading cache from file
my $cache = Bio::EnsEMBL::IdMapping::Cache->new(
  -LOGGER       => $logger,
  -CONF         => $conf,
);

foreach my $dbtype (qw(source target)) {
  foreach my $slice_name (@{ $cache->slice_names($dbtype) }) {
    $logger->log("\n");
    $cache->read_from_file('exons_by_id', "$dbtype.$slice_name");
  }
}

$cache->merge('exons_by_id');


# run in requested mode
my $mode = $conf->param('mode') || 'normal';
no strict 'refs';
&run_$mode;


# finish logfile
$logger->finish_log;


### END main ###


sub run_normal {
  # build scores
  &build_scores;
}


sub build_scores {
  my $esb = Bio::EnsEMBL::IdMapping::ExonScoreBuilder->new(
    -LOGGER       => $logger,
    -CONF         => $conf,
    -CACHE        => $cache
  );
  #my $tsb = Bio::EnsEMBL::IdMapping::TranscriptScoreBuilder->new(
  #my $gsb = Bio::EnsEMBL::IdMapping::GeneScoreBuilder->new(

  $exon_scores = $esb->score_exons;
  #$transcript_scores = $tsb->score_transcripts;
  #$gene_scores = $gsb->score_genes;
}


#
# test memory consumption of cache after merging. used for debugging.
#
sub log_cache_stats {
  $logger->log("\nCache memory usage:\n\n");

  my $s;
  my %keys;

  $keys{'cache'} = size($cache->{'cache'});

  foreach my $name (keys %{ $cache->{'cache'} }) {
    $keys{$name} = size($cache->{'cache'}->{$name});
    foreach my $type (keys %{ $cache->{'cache'}->{$name} }) {
      $keys{$type} = size($cache->{'cache'}->{$name}->{$type});
      $s += size($cache->{'cache'}->{$name}->{$type});
    }
  }

  my $ts = total_size($cache->{'cache'});

  my $fmt = "%-50s%12.0f\n";

  foreach my $k (sort { $keys{$a} <=> $keys{$b} } keys %keys) {
    $logger->log(sprintf($fmt, $k, $keys{$k}), 1);
  }
  $logger->log(sprintf($fmt, "total overhead", $s), 1);
  $logger->log(sprintf($fmt, "data", ($ts-$s)), 1);
  $logger->log(sprintf($fmt, "total", $ts)."\n", 1);

  # test
  my $i = 0;
  foreach my $eid (keys %{ $cache->get_by_name('exons_by_id', 'target') }) {
    last if ($i++ > 0);
    
    $logger->log("\nData object memory usage:\n\n");
    
    my $exon = $cache->get_by_key('exons_by_id', 'target', $eid);
    my $s1 = size($exon);
    my $ts1 = total_size($exon);

    $logger->log(sprintf($fmt, "object", $s1), 1);
    $logger->log(sprintf($fmt, "data", ($ts1-$s1)), 1);
    $logger->log(sprintf($fmt, "total", $ts1)."\n", 1);

    print $exon->stable_id."\n";
    #warn Data::Dumper::Dumper($exon);
  }
}


