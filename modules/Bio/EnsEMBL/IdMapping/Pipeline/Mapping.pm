=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2023] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

package Bio::EnsEMBL::IdMapping::Pipeline::Mapping;

use strict;
use warnings;
use Data::Dumper;
use File::Basename;

use base qw/Bio::EnsEMBL::Hive::Process/;
use Bio::EnsEMBL::Utils::ConfParser;
use Bio::EnsEMBL::Utils::Logger;
use Bio::EnsEMBL::IdMapping::Cache;
use Bio::EnsEMBL::Utils::ScriptUtils qw(path_append);
use Bio::EnsEMBL::IdMapping::StableIdMapper;
use Bio::EnsEMBL::IdMapping::ExonScoreBuilder;
use Bio::EnsEMBL::IdMapping::TranscriptScoreBuilder;
use Bio::EnsEMBL::IdMapping::GeneScoreBuilder;
use Bio::EnsEMBL::IdMapping::InternalIdMapper;
use Bio::EnsEMBL::IdMapping::Archiver;

my ($conf, $logger);
my ($cache, $stable_id_mapper);
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
my %mapping_types = ();

sub run {
  my ($self) = @_;

  my $config     = $self->param_required('config');
  my $mode       = $self->param_required('mode');
  my $logautoid  = $self->param('logautoid');
  my $logauto    = $self->param_required('logauto');

  # Parse configuration and commandline arguments
  $conf = new Bio::EnsEMBL::Utils::ConfParser(
    -SERVERROOT => dirname($config)."/../../..",
    -DEFAULT_CONF => $config
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

  $conf->param('mode', $mode);
  $conf->param('logauto', $logauto);

  # Set default logpath
  unless ($conf->param('logpath')) {
    $conf->param('logpath', path_append($conf->param('basedir'), 'log'));
  }
  
  $logger = new Bio::EnsEMBL::Utils::Logger(
    -LOGFILE      => $conf->param('logfile'),
    -LOGAUTO      => $conf->param('logauto'),
    -LOGAUTOBASE  => 'id_mapping',
    -LOGAUTOID    => $logautoid,
    -LOGPATH      => $conf->param('logpath'),
    -LOGAPPEND    => $conf->param('logappend'),
    -LOGLEVEL     => $conf->param('loglevel'),
    -IS_COMPONENT => 1,
  );

  # Initialise log
  $logger->init_log($conf->list_param_values);

  # Loading cache from file
  $cache = Bio::EnsEMBL::IdMapping::Cache->new(
    -LOGGER         => $logger,
    -CONF           => $conf,
    -LOAD_INSTANCE  => 1,
  );

  # Get a stable ID mapper
  $stable_id_mapper = Bio::EnsEMBL::IdMapping::StableIdMapper->new(
    -LOGGER       => $logger,
    -CONF         => $conf,
    -CACHE        => $cache
  );

  # Find out which entities we want to map
  foreach my $type ($conf->param('mapping_types')) {
    $mapping_types{$type} = 1;
  }

  if ($mode eq 'normal' || $mode eq 'mapping') {
    run_normal();
  } else {
    $logger->error("Invalid mode $mode.\n");
    $self->log()->error("Invalid mode $mode.");
  }

  # Finish logfile
  $logger->finish_log;
}

sub run_normal {
  my $self = shift;

  # Build scores
  build_scores();

  # Map stable IDs
  mapping();

  # Assign stable IDs and make creation and deletion events
  assign_stable_ids();

  # Generate similarity events
  $self->generate_similarity_events();

  # Dump existing stable_id_event table to file
  $self->dump_existing_events();

  # Create gene and peptide archive
  $self->archive($stable_id_mapper->mapping_session_id);

  # Upload table data files into db
  # upload_mapping_session_and_events();
  # upload_stable_ids();
  # upload_archive();

  # Final stats and mapping summary
  # analyse_results();
}

sub build_scores {
  # Get new ScoreBuilders for exons, transcripts and genes
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

  # Exon scoring
  $exon_scores = $esb->score_exons;
  
  # Transcript scoring
  $transcript_scores = $tsb->score_transcripts($exon_scores);
  
  # Gene scoring
  $gene_scores = $gsb->score_genes($transcript_scores);
}

sub mapping {
  # Get an internal ID mapper
  my $internal_id_mapper = Bio::EnsEMBL::IdMapping::InternalIdMapper->new(
    -LOGGER       => $logger,
    -CONF         => $conf,
    -CACHE        => $cache
  );

  # Map genes
  $gene_mappings = $internal_id_mapper->map_genes($gene_scores, $transcript_scores, $gsb);

  # Map transcripts
  if ($mapping_types{'transcript'} || $mapping_types{'exon'} || $mapping_types{'translation'}) {
    $transcript_mappings = $internal_id_mapper->map_transcripts($transcript_scores, $gene_mappings, $tsb);
  }

  # Map exons
  if ($mapping_types{'exon'}) {
    $exon_mappings = $internal_id_mapper->map_exons($exon_scores, $transcript_mappings, $esb);
  }

  # Map translations
  if ($mapping_types{'translation'}) {
    $translation_mappings = $internal_id_mapper->map_translations($transcript_mappings);
  }
}

sub assign_stable_ids {
  # Exons
  if ($mapping_types{'exon'}) {
    $stable_id_mapper->map_stable_ids($exon_mappings, 'exon');
  }

  # Transcripts
  if ($mapping_types{'transcript'}) {
    $stable_id_mapper->map_stable_ids($transcript_mappings, 'transcript');
  }

  # Translations
  if ($mapping_types{'translation'}) {
    $stable_id_mapper->map_stable_ids($translation_mappings, 'translation');
  }

  # Genes
  if ($mapping_types{'gene'}) {
    $stable_id_mapper->map_stable_ids($gene_mappings, 'gene');
  }


  # Dump mappings to file for debug purposes
  $stable_id_mapper->dump_debug_mappings;

  # Write stable_id_events to file
  $stable_id_mapper->write_stable_id_events('new');
}

sub generate_similarity_events {
  my $self = shift;

  $logger->info("Generating similarity events...\n", 0, 'stamped');
  $self->log()->info("Generating similarity events...");

  # Genes
  if ($mapping_types{'gene'}) {
    $logger->debug("genes\n", 1);
    $self->log()->debug("genes");
    $stable_id_mapper->generate_similarity_events($gene_mappings, $gene_scores, 'gene');
  }

  # Transcripts
  my $filtered_transcript_scores;
  if ($mapping_types{'transcript'} or $mapping_types{'translation'}) {
    $filtered_transcript_scores = $stable_id_mapper->filter_same_gene_transcript_similarities($transcript_scores);
  }

  if ($mapping_types{'transcript'}) {
    $logger->debug("transcripts\n", 1);
    $self->log()->debug("transcripts");
    $stable_id_mapper->generate_similarity_events($transcript_mappings, $filtered_transcript_scores, 'transcript');
  }

  # Translations
  if ($mapping_types{'translation'}) {
    $logger->debug("translations\n", 1);
    $self->log()->debug("translations");
    $stable_id_mapper->generate_translation_similarity_events($translation_mappings, $filtered_transcript_scores);
  }

  # Write stable_id_events to file
  $stable_id_mapper->write_stable_id_events('similarity');

  # write_retrofit_stable_id_events?? [todo]

  $logger->info("Done.\n\n", 0, 'stamped');
  $self->log()->info("Done");
}

sub dump_existing_events {
  my $self = shift;

  $logger->info("Dumping existing stable_id_events...\n", 0, 'stamped');
  $self->log()->info("Dumping existing stable_id_events...");
  
  my $i = $stable_id_mapper->dump_table_to_file('source', 'stable_id_event', 'stable_id_event_existing.txt', 1);

  $logger->info("Done writing $i entries.\n\n", 0, 'stamped');
  $self->log()->info("Done writing $i entries.");
}

sub archive {
  my $self = shift;
  my $mapping_session_id = shift;
  
  $logger->info("Create gene and peptide archive...\n", 0, 'stamped');
  $self->log()->info("Create gene and peptide archive...");
  
  # Get an Archiver
  my $archiver = Bio::EnsEMBL::IdMapping::Archiver->new(
    -LOGGER       => $logger,
    -CONF         => $conf,
    -CACHE        => $cache
  );

  # Create gene and peptide archive
  $archiver->create_archive($mapping_session_id);

  $logger->info("Done.\n\n", 0, 'stamped');
  $self->log()->info("Done.");

  # Dump existing archive tables to file
  $logger->info("Dumping existing gene and peptide archive...\n", 0, 'stamped');
  $self->log()->info("Dumping existing gene and peptide archive...");

  my $i = $archiver->dump_table_to_file('source', 'gene_archive', 'gene_archive_existing.txt', 1);
  my $j = $archiver->dump_table_to_file('source', 'peptide_archive', 'peptide_archive_existing.txt', 1);

  $logger->info("Done writing $i gene_archive and $j peptide_archive entries.\n\n", 0, 'stamped');
  $self->log()->info("Done writing $i gene_archive and $j peptide_archive entries.");
}

1;
