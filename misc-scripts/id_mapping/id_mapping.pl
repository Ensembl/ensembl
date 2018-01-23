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
use Bio::EnsEMBL::Utils::ScriptUtils qw(path_append);
use Bio::EnsEMBL::IdMapping::Cache;
use Bio::EnsEMBL::IdMapping::ExonScoreBuilder;
use Bio::EnsEMBL::IdMapping::TranscriptScoreBuilder;
use Bio::EnsEMBL::IdMapping::GeneScoreBuilder;
use Bio::EnsEMBL::IdMapping::InternalIdMapper;
use Bio::EnsEMBL::IdMapping::StableIdMapper;
use Bio::EnsEMBL::IdMapping::Archiver;
use Bio::EnsEMBL::IdMapping::ResultAnalyser;

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
  'basedir|basedir=s' => 1,
  'chromosomes|chr=s@' => 0,
  'region=s' => 0,
  'biotypes=s@' => 0,
  'biotypes_include=s@' => 0,
  'biotypes_exclude=s@' => 0,
  'min_exon_length|minexonlength=i' => 0,
  'exonerate_path|exoneratepath=s' => 1,
  'exonerate_threshold|exoneratethreshold=f' => 0,
  'exonerate_jobs|exoneratejobs=i' => 0,
  'exonerate_bytes_per_job|exoneratebytesperjob=f' => 0,
  'exonerate_extra_params|exonerateextraparams=s' => 0,
  'plugin_internal_id_mappers_gene=s@' => 0,
  'plugin_internal_id_mappers_transcript=s@' => 0,
  'plugin_internal_id_mappers_exon=s@' => 0,
  'mapping_types=s@' => 1,
  'plugin_stable_id_generator=s' => 0,
  'upload_events|uploadevents=s' => 0,
  'upload_stable_ids|uploadstableids=s' => 0,
  'upload_archive|uploadarchive=s' => 0,
);

# set default logpath
unless ($conf->param('logpath')) {
  $conf->param('logpath', path_append($conf->param('basedir'), 'log'));
}

# get log filehandle and print heading and parameters to logfile
my $logger = new Bio::EnsEMBL::Utils::Logger(
  -LOGFILE      => $conf->param('logfile'),
  -LOGAUTO      => $conf->param('logauto'),
  -LOGAUTOBASE  => 'id_mapping',
  -LOGAUTOID    => $conf->param('logautoid'),
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
my $exon_mappings;
my $transcript_mappings;
my $gene_mappings;
my $translation_mappings;

# loading cache from file
my $cache = Bio::EnsEMBL::IdMapping::Cache->new(
  -LOGGER         => $logger,
  -CONF           => $conf,
  -LOAD_INSTANCE  => 1,
);


# get a stable ID mapper
my $stable_id_mapper = Bio::EnsEMBL::IdMapping::StableIdMapper->new(
  -LOGGER       => $logger,
  -CONF         => $conf,
  -CACHE        => $cache
);


# find out which entities we want to map
my %mapping_types = ();
foreach my $type ($conf->param('mapping_types')) {
  $mapping_types{$type} = 1;
}


# run in requested mode
my $mode = $conf->param('mode') || 'normal';
if ( $mode eq 'mapping' ) { $mode = 'normal' }
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

  # assign stable IDs and make creation and deletion events
  &assign_stable_ids;

  # generate similarity events
  &generate_similarity_events;

  # dump existing stable_id_event table to file
  &dump_existing_events;

  # create gene and peptide archive
  &archive($stable_id_mapper->mapping_session_id);

  # upload table data files into db
  &upload_mapping_session_and_events;
  &upload_stable_ids;
  &upload_archive;

  # final stats and mapping summary
  &analyse_results;
}


sub run_upload {
  # upload table data files into db
  &upload_mapping_session_and_events;
  &upload_stable_ids;
  &upload_archive;
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
  
  # get an internal ID mapper
  my $internal_id_mapper = Bio::EnsEMBL::IdMapping::InternalIdMapper->new(
    -LOGGER       => $logger,
    -CONF         => $conf,
    -CACHE        => $cache
  );

  # map genes
  $gene_mappings = $internal_id_mapper->map_genes($gene_scores,
    $transcript_scores, $gsb);

  # map transcripts
  if ($mapping_types{'transcript'} or $mapping_types{'exon'} or
      $mapping_types{'translation'}) {
    $transcript_mappings = $internal_id_mapper->map_transcripts(
      $transcript_scores, $gene_mappings, $tsb);
  }

  # map exons
  if ($mapping_types{'exon'}) {
    $exon_mappings = $internal_id_mapper->map_exons($exon_scores,
      $transcript_mappings, $esb);
  }

  # map translations
  if ($mapping_types{'translation'}) {
    $translation_mappings = $internal_id_mapper->map_translations(
      $transcript_mappings);
  }
}


sub assign_stable_ids {
  
  #
  # assign stable IDs
  #

  # exons
  if ($mapping_types{'exon'}) {
    $stable_id_mapper->map_stable_ids($exon_mappings, 'exon');
  }

  # transcripts
  if ($mapping_types{'transcript'}) {
    $stable_id_mapper->map_stable_ids($transcript_mappings, 'transcript');
  }

  # translations
  if ($mapping_types{'translation'}) {
    $stable_id_mapper->map_stable_ids($translation_mappings, 'translation');
  }

  # genes
  if ($mapping_types{'gene'}) {
    $stable_id_mapper->map_stable_ids($gene_mappings, 'gene');
  }


  # dump mappings to file for debug purposes
  $stable_id_mapper->dump_debug_mappings;

  # write stable_id_events to file
  $stable_id_mapper->write_stable_id_events('new');

}


sub generate_similarity_events {
  
  $logger->info("Generating similarity events...\n", 0, 'stamped');

  # genes
  if ($mapping_types{'gene'}) {
    $logger->debug("genes\n", 1);
    $stable_id_mapper->generate_similarity_events($gene_mappings, $gene_scores,
      'gene');
  }

  # transcripts
  my $filtered_transcript_scores;
  if ($mapping_types{'transcript'} or $mapping_types{'translation'}) {
    $filtered_transcript_scores =
      $stable_id_mapper->filter_same_gene_transcript_similarities(
        $transcript_scores);
  }

  if ($mapping_types{'transcript'}) {
    $logger->debug("transcripts\n", 1);
    $stable_id_mapper->generate_similarity_events($transcript_mappings,
      $filtered_transcript_scores, 'transcript');
  }

  # translations
  if ($mapping_types{'translation'}) {
    $logger->debug("translations\n", 1);
    $stable_id_mapper->generate_translation_similarity_events(
      $translation_mappings, $filtered_transcript_scores);
  }

  # write stable_id_events to file
  $stable_id_mapper->write_stable_id_events('similarity');

  # write_retrofit_stable_id_events?? [todo]

  $logger->info("Done.\n\n", 0, 'stamped');
}


sub dump_existing_events {
  $logger->info("Dumping existing stable_id_events...\n", 0, 'stamped');
  
  my $i = $stable_id_mapper->dump_table_to_file('source', 'stable_id_event',
    'stable_id_event_existing.txt', 1);

  $logger->info("Done writing $i entries.\n\n", 0, 'stamped');
}


sub archive {
  my $mapping_session_id = shift;
  
  $logger->info("Create gene and peptide archive...\n", 0, 'stamped');
  
  # get an Archiver
  my $archiver = Bio::EnsEMBL::IdMapping::Archiver->new(
    -LOGGER       => $logger,
    -CONF         => $conf,
    -CACHE        => $cache
  );

  # create gene and peptide archive
  $archiver->create_archive($mapping_session_id);

  $logger->info("Done.\n\n", 0, 'stamped');

  # dump existing archive tables to file
  $logger->info("Dumping existing gene and peptide archive...\n", 0, 'stamped');

  my $i = $archiver->dump_table_to_file('source', 'gene_archive',
    'gene_archive_existing.txt', 1);
  my $j = $archiver->dump_table_to_file('source', 'peptide_archive',
    'peptide_archive_existing.txt', 1);

  $logger->info("Done writing $i gene_archive and $j peptide_archive entries.\n\n", 0, 'stamped');
}


sub upload_mapping_session_and_events {
  if ($conf->is_true('upload_events') and ! $conf->param('dry_run')) {
    
    $logger->info("Uploading mapping_session and stable_id_event tables...\n");

    my $i = 0;
    my $j = 0;
    
    $logger->info("mapping_session...\n", 1);
    $i += $stable_id_mapper->upload_file_into_table('target', 'mapping_session',
      'mapping_session.txt');
    $logger->info("$i\n", 1);
    
    $logger->info("stable_id_event...\n", 1);
    $j += $stable_id_mapper->upload_file_into_table('target', 'stable_id_event',
      'stable_id_event_existing.txt');
    $j += $stable_id_mapper->upload_file_into_table('target', 'stable_id_event',
      'stable_id_event_new.txt', 1);
    $j += $stable_id_mapper->upload_file_into_table('target', 'stable_id_event',
      'stable_id_event_similarity.txt', 1);
    $logger->info("$j\n", 1);
    
    $logger->info("Done.\n\n");
  
  } else {
    $logger->info("Stable ID event and mapping session tables not uploaded.\n\n");
  }
}


sub upload_stable_ids {
  if ($conf->is_true('upload_stable_ids') and ! $conf->param('dry_run')) {
    
    $logger->info("Uploading stable ID tables...\n");
    
    foreach my $t ($conf->param('mapping_types')) {
      $logger->info("${t}_stable_id...\n", 1);
      my $i = $stable_id_mapper->upload_file_into_table('target',
        "${t}_stable_id", "${t}_stable_id.txt");
      $logger->info("$i\n", 1);
    }
    
    $logger->info("Done.\n\n");
  
  } else {
    $logger->info("Stable ID tables not uploaded.\n\n");
  }
}


sub upload_archive {
  if ($conf->is_true('upload_archive') and ! $conf->param('dry_run')) {
    
    $logger->info("Uploading gene and peptide tables...\n");
    
    foreach my $t (qw(gene peptide)) {
      $logger->info("${t}_archive...\n", 1);
      my $i = 0;
      $i += $stable_id_mapper->upload_file_into_table('target', "${t}_archive",
        "${t}_archive_existing.txt", 1);
      $i += $stable_id_mapper->upload_file_into_table('target', "${t}_archive",
        "${t}_archive_new.txt", 1);
      $logger->info("$i\n", 1);
    }
    
    $logger->info("Done.\n\n");
  
  } else {
    $logger->info("Gene and peptide archive tables not uploaded.\n\n");
  }
}


sub analyse_results {
    
  $logger->info("Analysing results...\n", 0, 'stamped');
  
  # get a result analyser
  my $analyser = Bio::EnsEMBL::IdMapping::ResultAnalyser->new(
    -LOGGER       => $logger,
    -CONF         => $conf,
    -CACHE        => $cache
  );

  # analyse results
  $analyser->analyse($gene_mappings,
    $stable_id_mapper->get_all_stable_id_events('similarity'));

  # write results to file
  $analyser->write_results_to_file;

  # create click lists
  $analyser->create_clicklist;
  
  $logger->info("Done.\n\n", 0, 'stamped');

  # mapping summary
  $logger->info("Creating mapping summary...\n", 0, 'stamped');
  $analyser->create_mapping_summary;
  $logger->info("Done.\n", 0, 'stamped');
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


