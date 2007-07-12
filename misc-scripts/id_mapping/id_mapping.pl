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
use Bio::EnsEMBL::IdMapping::Cache;
use Bio::EnsEMBL::IdMapping::ExonScoreBuilder;
use Bio::EnsEMBL::IdMapping::TranscriptScoreBuilder;
use Bio::EnsEMBL::IdMapping::GeneScoreBuilder;
use Bio::EnsEMBL::IdMapping::InternalIdMapper;

#use Devel::Size qw(size total_size);
#use Data::Dumper;
#$Data::Dumper::Indent = 1;

# parse configuration and commandline arguments
my $conf = new Bio::EnsEMBL::Utils::ConfParser(
  -SERVERROOT => "$Bin/../../..",
  -DEFAULT_CONF => "$Bin/default.conf"
);

$conf->parse_options(
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
);

# get log filehandle and print heading and parameters to logfile
my $logger = new Bio::EnsEMBL::Utils::Logger(
  -LOGFILE      => $conf->param('logfile'),
  -LOGPATH      => $conf->param('logpath'),
  -LOGAPPEND    => $conf->param('logappend'),
  -LOGLEVEL     => $conf->param('loglevel'),
  -IS_COMPONENT => $conf->param('is_component'),
);

# initialise log
$logger->init_log($conf->list_param_values);


# instance variables
my $esb;
my $tsb;
my $gsb;
my $exon_scores;
my $transcript_scores;
my $gene_scores;


# loading cache from file
my $cache = Bio::EnsEMBL::IdMapping::Cache->new(
  -LOGGER       => $logger,
  -CONF         => $conf,
);
$cache->read_instance_from_file;


# run in requested mode
my $mode = $conf->param('mode') || 'normal';
my $run = "run_$mode";
no strict 'refs';
&$run;


# finish logfile
$logger->finish_log;


### END main ###


sub run_normal {
  # build scores
  &build_scores;

  # map stable IDs
  &map;
}


sub build_scores {
  # get new ScoreBuilders for exons, transcripts and genes
  $esb = Bio::EnsEMBL::IdMapping::ExonScoreBuilder->new(
    -LOGGER       => $logger,
    -CONF         => $conf,
    -CACHE        => $cache
  );
  $tsb = Bio::EnsEMBL::IdMapping::TranscriptScoreBuilder->new(
    -LOGGER       => $logger,
    -CONF         => $conf,
    -CACHE        => $cache
  );
  $gsb = Bio::EnsEMBL::IdMapping::GeneScoreBuilder->new(
    -LOGGER       => $logger,
    -CONF         => $conf,
    -CACHE        => $cache
  );

  # exon scoring
  $exon_scores = $esb->score_exons;
  
  # transcript scoring
  $transcript_scores = $tsb->score_transcripts($exon_scores);
  
  # gene scoring
  $gene_scores = $gsb->score_genes($transcript_scores);
}


sub map {
  
  # get a mapper
  my $mapper = Bio::EnsEMBL::IdMapping::InternalIdMapper->new(
    -LOGGER       => $logger,
    -CONF         => $conf,
    -CACHE        => $cache
  );

  # map genes
  $mapper->map_genes($gene_scores, $transcript_scores, $gsb);

  # map transcripts

  # map exons

  # map translations
}


#
# test memory consumption of cache after merging. used for debugging.
#
sub log_cache_stats {
  $logger->info("\nCache memory usage:\n\n");

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
    $logger->info(sprintf($fmt, $k, $keys{$k}), 1);
  }
  $logger->info(sprintf($fmt, "total overhead", $s), 1);
  $logger->info(sprintf($fmt, "data", ($ts-$s)), 1);
  $logger->info(sprintf($fmt, "total", $ts)."\n", 1);

  # test
  my $i = 0;
  foreach my $eid (keys %{ $cache->get_by_name('exons_by_id', 'target') }) {
    last if ($i++ > 0);
    
    $logger->info("\nData object memory usage:\n\n");
    
    my $exon = $cache->get_by_key('exons_by_id', 'target', $eid);
    my $s1 = size($exon);
    my $ts1 = total_size($exon);

    $logger->info(sprintf($fmt, "object", $s1), 1);
    $logger->info(sprintf($fmt, "data", ($ts1-$s1)), 1);
    $logger->info(sprintf($fmt, "total", $ts1)."\n", 1);

    print $exon->stable_id."\n";
    #warn Data::Dumper::Dumper($exon);
  }
}


