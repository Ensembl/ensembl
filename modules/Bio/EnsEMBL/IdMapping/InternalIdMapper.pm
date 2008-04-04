package Bio::EnsEMBL::IdMapping::InternalIdMapper;

=head1 NAME


=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS


=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http:#www.ensembl.org/info/about/code_licence.html for details.

=head1 AUTHOR

Patrick Meidl <meidl@ebi.ac.uk>, Ensembl core API team

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<ensembl-dev@ebi.ac.uk>

=cut


use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::IdMapping::BaseObject;
our @ISA = qw(Bio::EnsEMBL::IdMapping::BaseObject);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::ScriptUtils qw(path_append);
use Bio::EnsEMBL::IdMapping::Entry;
use Bio::EnsEMBL::IdMapping::MappingList;
use Bio::EnsEMBL::IdMapping::SyntenyFramework;


# scores are considered the same if (2.0 * (s1-s2))/(s1 + s2) < this
use constant SIMILAR_SCORE_RATIO => 0.01;

    
sub map_genes {
  my $self = shift;
  my $gene_scores = shift;
  my $transcript_scores = shift;
  my $gsb = shift;

  # argument checks
  unless ($gene_scores and
          $gene_scores->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix')) {
    throw('Need a gene Bio::EnsEMBL::IdMapping::ScoredMappingMatrix.');
  }
  
  unless ($transcript_scores and
          $transcript_scores->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix')) {
    throw('Need a transcript Bio::EnsEMBL::IdMapping::ScoredMappingMatrix.');
  }

  unless ($gsb and
          $gsb->isa('Bio::EnsEMBL::IdMapping::GeneScoreBuilder')) {
    throw('Need a Bio::EnsEMBL::IdMapping::GeneScoreBuilder.');
  }
  
  $self->logger->info("== Internal ID mapping for genes...\n\n", 0, 'stamped');

  my $dump_path = path_append($self->conf->param('dumppath'), 'mapping');

  my $mappings = Bio::EnsEMBL::IdMapping::MappingList->new(
    -DUMP_PATH   => $dump_path,
    -CACHE_FILE  => 'gene_mappings.ser',
  );

  my $mapping_cache = $mappings->cache_file;

  if (-s $mapping_cache) {
    
    # read from file
    $self->logger->info("Reading gene mappings from file...\n", 0, 'stamped');
    $self->logger->debug("Cache file $mapping_cache.\n", 1);
    $mappings->read_from_file;
    $self->logger->info("Done.\n\n", 0, 'stamped');
    
  } else {
    
    # create gene mappings
    $self->logger->info("No gene mappings found. Will calculate them now.\n");

    #
    # basic mapping
    #
    $self->logger->info("Basic gene mapping...\n", 0, 'stamped');

    my $mappings0 = $self->basic_mapping($gene_scores, 'gene_mappings0');

    my $gene_scores1 = $gsb->create_shrinked_matrix($gene_scores, $mappings0,
      'gene_matrix1');


    #
    # build the synteny from unambiguous mappings
    #
    unless ($gene_scores1->loaded) {
      $self->logger->info("Synteny Framework building...\n", 0, 'stamped');
      my $sf = Bio::EnsEMBL::IdMapping::SyntenyFramework->new(
        -DUMP_PATH    => $dump_path,
        -CACHE_FILE   => 'synteny_framework.ser',
        -LOGGER       => $self->logger,
        -CONF         => $self->conf,
        -CACHE        => $self->cache,
      );
      $sf->build_synteny($mappings0);

      # use it to rescore the genes
      $self->logger->info("\nSynteny assisted mapping...\n", 0, 'stamped');
      $gene_scores1 = $sf->rescore_gene_matrix_lsf($gene_scores1);

      # checkpoint
      $gene_scores1->write_to_file;
    }

    my $mappings1 = $self->basic_mapping($gene_scores1, 'gene_mappings1');
    
    my $gene_scores2 = $gsb->create_shrinked_matrix($gene_scores1, $mappings1,
      'gene_matrix2');
    

    #
    # rescore with simple scoring function and try again
    #
    $self->logger->info("Retry with simple best transcript score...\n", 0, 'stamped');
    
    unless ($gene_scores2->loaded) {
      $gsb->simple_gene_rescore($gene_scores2, $transcript_scores);
      $gene_scores2->write_to_file;
    }
    
    my $mappings2 = $self->basic_mapping($gene_scores2, 'gene_mappings2');
    
    my $gene_scores3 = $gsb->create_shrinked_matrix($gene_scores2, $mappings2,
      'gene_matrix3');


    #
    # rescore by penalising scores between genes with different biotypes  
    #
    $self->logger->info("Retry with biotype disambiguation...\n", 0, 'stamped');
    
    unless ($gene_scores3->loaded) {
      $gsb->biotype_gene_rescore($gene_scores3);
      $gene_scores3->write_to_file;
    }

    my $mappings3 = $self->basic_mapping($gene_scores3, 'gene_mappings3');
    
    my $gene_scores4 = $gsb->create_shrinked_matrix($gene_scores3, $mappings3,
      'gene_matrix4');


    #
    # selectively rescore by penalising scores between genes with different
    # internalIDs  
    #
    $self->logger->info("Retry with internalID disambiguation...\n", 0, 'stamped');
    
    unless ($gene_scores4->loaded) {
      $gsb->internal_id_rescore($gene_scores4);
      $gene_scores4->write_to_file;
    }

    my $mappings4 = $self->basic_mapping($gene_scores4, 'gene_mappings4');
    
    my $remaining_gene_scores = $gsb->create_shrinked_matrix(
      $gene_scores4, $mappings4, 'remaining_gene_matrix');


    #
    # report remaining ambiguities
    #
    $self->logger->info($remaining_gene_scores->get_source_count.
      " source genes are ambiguous with ".
      $remaining_gene_scores->get_target_count." target genes.\n\n");

    $self->log_ambiguous($remaining_gene_scores, 'gene');

    
    #
    # merge mappings and write to file
    #
    $mappings->add_all($mappings0, $mappings1, $mappings2, $mappings3,
                       $mappings4);

    $mappings->write_to_file;

    if ($self->logger->loglevel eq 'debug') {
      $mappings->log('gene', $self->conf->param('dumppath'));
    }

    $self->logger->info("Done.\n\n", 0, 'stamped');

  }

  return $mappings;
}


sub map_transcripts {
  my $self = shift;
  my $transcript_scores = shift;
  my $gene_mappings = shift;
  my $tsb = shift;

  # argument checks
  unless ($transcript_scores and
      $transcript_scores->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix')) {
    throw('Need a transcript Bio::EnsEMBL::IdMapping::ScoredMappingMatrix.');
  }

  unless ($gene_mappings and
          $gene_mappings->isa('Bio::EnsEMBL::IdMapping::MappingList')) {
    throw('Need a gene Bio::EnsEMBL::IdMapping::MappingList.');
  }
  
  unless ($tsb and
          $tsb->isa('Bio::EnsEMBL::IdMapping::TranscriptScoreBuilder')) {
    throw('Need a Bio::EnsEMBL::IdMapping::TranscriptScoreBuilder.');
  }
  
  $self->logger->info("== Internal ID mapping for transcripts...\n\n", 0, 'stamped');

  my $dump_path = path_append($self->conf->param('dumppath'), 'mapping');

  my $mappings = Bio::EnsEMBL::IdMapping::MappingList->new(
    -DUMP_PATH   => $dump_path,
    -CACHE_FILE  => 'transcript_mappings.ser',
  );

  my $mapping_cache = $mappings->cache_file;

  if (-s $mapping_cache) {
    
    # read from file
    $self->logger->info("Reading transcript mappings from file...\n", 0,
      'stamped');
    $self->logger->debug("Cache file $mapping_cache.\n", 1);
    $mappings->read_from_file;
    $self->logger->info("Done.\n\n", 0, 'stamped');
    
  } else {
    
    # create transcript mappings
    $self->logger->info("No transcript mappings found. Will calculate them now.\n");

    #
    # basic mapping
    #
    $self->logger->info("Basic transcript mapping...\n", 0, 'stamped');

    my $mappings0 = $self->basic_mapping($transcript_scores,
      'transcript_mappings0');

    my $transcript_scores1 = $tsb->create_shrinked_matrix(
      $transcript_scores, $mappings0, 'transcript_matrix1');


    #
    # handle cases with exact match but different translation
    #
    $self->logger->info("Exact Transcript non-exact Translation...\n", 0, 'stamped');
    
    unless ($transcript_scores1->loaded) {
      $tsb->different_translation_rescore($transcript_scores1);
      $transcript_scores1->write_to_file;
    }
    
    my $mappings1 = $self->basic_mapping($transcript_scores1,
      'transcript_mappings1');
    
    my $transcript_scores2 = $tsb->create_shrinked_matrix(
      $transcript_scores1, $mappings1, 'transcript_matrix2');


    #
    # reduce score for mappings of transcripts which do not belong to mapped
    # genes
    #
    $self->logger->info("Transcripts in mapped genes...\n", 0, 'stamped');
    
    unless ($transcript_scores2->loaded) {
    $tsb->non_mapped_gene_rescore($transcript_scores2, $gene_mappings);
      $transcript_scores2->write_to_file;
    }
    
    my $mappings2 = $self->basic_mapping($transcript_scores2,
      'transcript_mappings2');
    
    my $transcript_scores3 = $tsb->create_shrinked_matrix(
      $transcript_scores2, $mappings2, 'transcript_matrix3');


    #
    # selectively rescore by penalising scores between transcripts with
    # different internalIDs  
    #
    $self->logger->info("Retry with internalID disambiguation...\n", 0, 'stamped');
    
    unless ($transcript_scores3->loaded) {
      $tsb->internal_id_rescore($transcript_scores3);
      $transcript_scores3->write_to_file;
    }

    my $mappings3 = $self->basic_mapping($transcript_scores3,
      'transcript_mappings3');
    
    my $transcript_scores4 = $tsb->create_shrinked_matrix(
      $transcript_scores3, $mappings3, 'transcript_matrix4');


    #
    # handle ambiguities between transcripts in single genes
    #
    $self->logger->info("Transcripts in single genes...\n", 0, 'stamped');
    
    unless ($transcript_scores4->loaded) {
      $transcript_scores4->write_to_file;
    }
    
    my $mappings4 = $self->same_gene_transcript_mapping($transcript_scores4,
      'transcript_mappings4');

    my $remaining_transcript_scores = $tsb->create_shrinked_matrix(
      $transcript_scores4, $mappings4, 'transcript_matrix5');


    #
    # report remaining ambiguities
    #
    $self->logger->info($remaining_transcript_scores->get_source_count.
      " source transcripts are ambiguous with ".
      $remaining_transcript_scores->get_target_count." target transcripts.\n\n");

    $self->log_ambiguous($remaining_transcript_scores, 'transcript');

    
    #
    # merge mappings and write to file
    #
    $mappings->add_all($mappings0, $mappings1, $mappings2, $mappings3,
                       $mappings4);

    $mappings->write_to_file;

    if ($self->logger->loglevel eq 'debug') {
      $mappings->log('transcript');
    }

    $self->logger->info("Done.\n\n", 0, 'stamped');

  }

  return $mappings;

}


sub map_exons {
  my $self = shift;
  my $exon_scores = shift;
  my $transcript_mappings = shift;
  my $esb = shift;

  # argument checks
  unless ($exon_scores and
      $exon_scores->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix')) {
    throw('Need a Bio::EnsEMBL::IdMapping::ScoredMappingMatrix of exons.');
  }

  unless ($transcript_mappings and
          $transcript_mappings->isa('Bio::EnsEMBL::IdMapping::MappingList')) {
    throw('Need a Bio::EnsEMBL::IdMapping::MappingList of transcripts.');
  }
  
  unless ($esb and
          $esb->isa('Bio::EnsEMBL::IdMapping::ExonScoreBuilder')) {
    throw('Need a Bio::EnsEMBL::IdMapping::ExonScoreBuilder.');
  }
  
  $self->logger->info("== Internal ID mapping for exons...\n\n", 0, 'stamped');

  my $dump_path = path_append($self->conf->param('dumppath'), 'mapping');

  my $mappings = Bio::EnsEMBL::IdMapping::MappingList->new(
    -DUMP_PATH   => $dump_path,
    -CACHE_FILE  => 'exon_mappings.ser',
  );

  my $mapping_cache = $mappings->cache_file;

  if (-s $mapping_cache) {
    
    # read from file
    $self->logger->info("Reading exon mappings from file...\n", 0,
      'stamped');
    $self->logger->debug("Cache file $mapping_cache.\n", 1);
    $mappings->read_from_file;
    $self->logger->info("Done.\n\n", 0, 'stamped');
    
  } else {
    
    # create exon mappings
    $self->logger->info("No exon mappings found. Will calculate them now.\n");

    #
    # basic mapping
    #
    $self->logger->info("Basic exon mapping...\n", 0, 'stamped');

    my $mappings0 = $self->basic_mapping($exon_scores, 'exon_mappings0');

    my $exon_scores1 = $esb->create_shrinked_matrix( $exon_scores, $mappings0,
      'exon_matrix1');
    

    #
    # reduce score for mappings of exons which do not belong to mapped
    # transcripts
    #
    $self->logger->info("Exons in mapped transcripts...\n", 0, 'stamped');
    
    unless ($exon_scores1->loaded) {
      $esb->non_mapped_transcript_rescore($exon_scores1, $transcript_mappings);
      $exon_scores1->write_to_file;
    }
    
    my $mappings1 = $self->basic_mapping($exon_scores1, 'exon_mappings1');
    
    my $remaining_exon_scores = $esb->create_shrinked_matrix(
      $exon_scores1, $mappings1, 'exon_matrix2');


    #
    # report remaining ambiguities
    #
    $self->logger->info($remaining_exon_scores->get_source_count.
      " source exons are ambiguous with ".
      $remaining_exon_scores->get_target_count." target exons.\n\n");

    $self->log_ambiguous($remaining_exon_scores, 'exon');

    
    #
    # merge mappings and write to file
    #
    $mappings->add_all($mappings0, $mappings1);

    $mappings->write_to_file;

    if ($self->logger->loglevel eq 'debug') {
      $mappings->log('exon');
    }

    $self->logger->info("Done.\n\n", 0, 'stamped');

  }

  return $mappings;

}


sub map_translations {
  my $self = shift;
  my $transcript_mappings = shift;

  # argument checks
  unless ($transcript_mappings and
          $transcript_mappings->isa('Bio::EnsEMBL::IdMapping::MappingList')) {
    throw('Need a Bio::EnsEMBL::IdMapping::MappingList of transcripts.');
  }
  
  $self->logger->info("== Internal ID mapping for translations...\n\n", 0, 'stamped');

  my $dump_path = path_append($self->conf->param('dumppath'), 'mapping');

  my $mappings = Bio::EnsEMBL::IdMapping::MappingList->new(
    -DUMP_PATH   => $dump_path,
    -CACHE_FILE  => 'translation_mappings.ser',
  );

  my $mapping_cache = $mappings->cache_file;

  if (-s $mapping_cache) {
    
    # read from file
    $self->logger->info("Reading translation mappings from file...\n", 0,
      'stamped');
    $self->logger->debug("Cache file $mapping_cache.\n", 1);
    $mappings->read_from_file;
    $self->logger->info("Done.\n\n", 0, 'stamped');
    
  } else {
    
    # create translation mappings
    $self->logger->info("No translation mappings found. Will calculate them now.\n");

    $self->logger->info("Translation mapping...\n", 0, 'stamped');

    #
    # map translations for mapped transcripts
    #
    my $i = 0;

    foreach my $entry (@{ $transcript_mappings->get_all_Entries }) {

      my $source_tl = $self->cache->get_by_key('transcripts_by_id',
        'source', $entry->source)->translation;
      my $target_tl = $self->cache->get_by_key('transcripts_by_id',
        'target', $entry->target)->translation;

      if ($source_tl and $target_tl) {
      
        # add mapping for the translations; note that the score is taken from
        # the transcript mapping
        my $tl_entry = Bio::EnsEMBL::IdMapping::Entry->new_fast([
          $source_tl->id, $target_tl->id, $entry->score
        ]);
        $mappings->add_Entry($tl_entry);
      
      } else {
        $i++;
      }

    }

    $self->logger->debug("Skipped transcripts without translation: $i\n", 1);
    $self->logger->info("New mappings: ".$mappings->get_entry_count."\n\n");

    $mappings->write_to_file;

    if ($self->logger->loglevel eq 'debug') {
      $mappings->log('translation');
    }

    $self->logger->info("Done.\n\n", 0, 'stamped');

  }

  return $mappings;

}


#
# find the highest unambiguous score for all sources and targets in a scoring
# matrix
#
sub basic_mapping {
  my $self = shift;
  my $matrix = shift;
  my $mapping_name = shift;

  # argument checks
  unless ($matrix and
          $matrix->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix')) {
    throw('Need a Bio::EnsEMBL::IdMapping::ScoredMappingMatrix.');
  }

  throw('Need a name for serialising the mapping.') unless ($mapping_name);

  # Create a new MappingList object. Specify AUTO_LOAD to load serialised
  # existing mappings if found
  my $dump_path = path_append($self->conf->param('dumppath'), 'mapping');
  
  my $mappings = Bio::EnsEMBL::IdMapping::MappingList->new(
    -DUMP_PATH   => $dump_path,
    -CACHE_FILE  => "${mapping_name}.ser",
    -AUTO_LOAD   => 1,
  );
  
  # checkpoint test: return a previously stored MappingList
  if ($mappings->loaded) {
    $self->logger->info("Read existing mappings from ${mapping_name}.ser.\n");
    return $mappings;
  }

  my $sources_done = {};
  my $targets_done = {};

  # sort scoring matrix entries by descending score
  my @sorted_entries = sort { $b->score <=> $a->score }
    @{ $matrix->get_all_Entries };

  # debug
  my $idx = substr($mapping_name, -1);

  while (my $entry = shift(@sorted_entries)) {
    
    #$self->logger->debug("\nxxx$idx ".$entry->to_string." ");
    
    # we already found a mapping for either source or target
    next if ($sources_done->{$entry->source} or
             $targets_done->{$entry->target});
    
    #$self->logger->debug('d');
    
    # there's a better mapping for either source or target
    next if ($self->higher_score_exists($entry, $matrix, $sources_done,
      $targets_done));
      
    #$self->logger->debug('h');

    # check for ambiguous mappings; they are dealt with later
    my $other_sources = [];
    my $other_targets = [];

    if ($self->ambiguous_mapping($entry, $matrix, $other_sources, $other_targets)) {
      #$self->logger->debug('a');
      
      $other_sources = $self->filter_sources($other_sources, $sources_done);
      $other_targets = $self->filter_targets($other_targets, $targets_done);

      next if (scalar(@$other_sources) or scalar(@$other_targets));
    }
    
    #$self->logger->debug('A');

    # this is the best mapping, add it
    $mappings->add_Entry($entry);

    $sources_done->{$entry->source} = 1;
    $targets_done->{$entry->target} = 1;
  }

  # create checkpoint
  $mappings->write_to_file;

  return $mappings;
}


sub higher_score_exists {
  my ($self, $entry, $matrix, $sources_done, $targets_done) = @_;

  my $source = $entry->source;
  my $target = $entry->target;
  my $score = $entry->score;

  foreach my $other_source (@{ $matrix->get_sources_for_target($target) }) {
    if ($other_source != $source and !$sources_done->{$other_source} and
        $score < $matrix->get_score($other_source, $target)) {
          return 1;
    }
  }

  foreach my $other_target (@{ $matrix->get_targets_for_source($source) }) {
    if ($other_target != $target and !$targets_done->{$other_target} and
        $score < $matrix->get_score($source, $other_target)) {
          return 1;
    }
  }

  return 0;
}


#
# find ambiguous mappings (see scores_similar() for definition)
#
sub ambiguous_mapping {
  my ($self, $entry, $matrix, $other_sources, $other_targets) = @_;

  my $source = $entry->source;
  my $target = $entry->target;
  my $score = $entry->score;

  my $retval = 0;

  foreach my $other_source (@{ $matrix->get_sources_for_target($target) }) {
    my $other_score = $matrix->get_score($other_source, $target);
    
    if ($other_source != $source and
      ($self->scores_similar($score, $other_score) or $score < $other_score)) {
        $retval = 1;
        push @{ $other_sources }, $other_source;
    }
  }

  foreach my $other_target (@{ $matrix->get_targets_for_source($source) }) {
    my $other_score = $matrix->get_score($source, $other_target);

    if ($other_target != $target and
      ($self->scores_similar($score, $other_score) or $score < $other_score)) {
        $retval = 1;
        push @{ $other_targets }, $other_target;
    }
  }

  return $retval;
}


# 
# rule for similarity taken from java code...
#
sub scores_similar {
  my ($self, $s1, $s2) = @_;

  # always give priority to exact matches over very similar ones
  return 0 if ($s1 == 1 and $s2 < 1);

  my $diff = $s1 -$s2;
  $diff = -$diff if ($diff < 0);
  
  my $pc = 2 * $diff / ($s1 + $s2);
  
  return ($pc < SIMILAR_SCORE_RATIO);
}


sub filter_sources {
  my ($self, $other_sources, $sources_done) = @_;

  unless (scalar(@$other_sources) and scalar(keys %$sources_done)) {
    return $other_sources;
  }

  my @tmp = ();

  foreach my $e (@{ $other_sources }) {
    push @tmp, $e unless ($sources_done->{$e}); 
  }

  return \@tmp;
}


sub filter_targets {
  my ($self, $other_targets, $targets_done) = @_;

  unless (scalar(@{ $other_targets }) and scalar(keys %$targets_done)) {
    return $other_targets;
  }

  my @tmp = ();

  foreach my $e (@{ $other_targets }) {
    push @tmp, $e unless ($targets_done->{$e}); 
  }

  return \@tmp;
}


#
# modified basic mapper that maps transcripts that are ambiguous within one gene
#
sub same_gene_transcript_mapping {
  my $self = shift;
  my $matrix = shift;
  my $mapping_name = shift;

  # argument checks
  unless ($matrix and
          $matrix->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix')) {
    throw('Need a Bio::EnsEMBL::IdMapping::ScoredMappingMatrix.');
  }

  throw('Need a name for serialising the mapping.') unless ($mapping_name);

  # Create a new MappingList object. Specify AUTO_LOAD to load serialised
  # existing mappings if found
  my $dump_path = path_append($self->conf->param('dumppath'), 'mapping');
  
  my $mappings = Bio::EnsEMBL::IdMapping::MappingList->new(
    -DUMP_PATH   => $dump_path,
    -CACHE_FILE  => "${mapping_name}.ser",
    -AUTO_LOAD   => 1,
  );
  
  # checkpoint test: return a previously stored MappingList
  if ($mappings->loaded) {
    $self->logger->info("Read existing mappings from ${mapping_name}.ser.\n");
    return $mappings;
  }

  my $sources_done = {};
  my $targets_done = {};

  # sort scoring matrix entries by descending score
  my @sorted_entries = sort { $b->score <=> $a->score ||
    $a->source <=> $b->source || $a->target <=> $b->target }
      @{ $matrix->get_all_Entries };

  while (my $entry = shift(@sorted_entries)) {
    
    # $self->logger->debug("\nxxx4 ".$entry->to_string." ");

    # we already found a mapping for either source or target yet
    next if ($sources_done->{$entry->source} or
             $targets_done->{$entry->target});

    #$self->logger->debug('d');

    my $other_sources = [];
    my $other_targets = [];
    my %source_genes = ();
    my %target_genes = ();

    if ($self->ambiguous_mapping($entry, $matrix, $other_sources, $other_targets)) {
      #$self->logger->debug('a');

      $other_sources = $self->filter_sources($other_sources, $sources_done);
      $other_targets = $self->filter_targets($other_targets, $targets_done);

      $source_genes{$self->cache->get_by_key('genes_by_transcript_id',
        'source', $entry->source)} = 1;
      $target_genes{$self->cache->get_by_key('genes_by_transcript_id',
        'target', $entry->target)} = 1;

      foreach my $other_source (@{ $other_sources }) {
        $source_genes{$self->cache->get_by_key('genes_by_transcript_id',
          'source', $other_source)} = 1;
      }
        
      foreach my $other_target (@{ $other_targets }) {
        $target_genes{$self->cache->get_by_key('genes_by_transcript_id',
          'target', $other_target)} = 1;
      }
      
      # only add mapping if only one source and target gene involved
      if (scalar(keys %source_genes) == 1 and scalar(keys %target_genes) == 1) {
        #$self->logger->debug('O');
        $mappings->add_Entry($entry);
      }

    } else {
      #$self->logger->debug('A');

      # this is the best mapping, add it
      $mappings->add_Entry($entry);
    }

    $sources_done->{$entry->source} = 1;
    $targets_done->{$entry->target} = 1;
  }

  # create checkpoint
  $mappings->write_to_file;

  return $mappings;
}


sub log_ambiguous {
  my $self = shift;
  my $matrix = shift;
  my $type = shift;

  unless ($matrix and
          $matrix->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix')) {
    throw('Need a Bio::EnsEMBL::IdMapping::ScoredMappingMatrix.');
  }
  
  # create dump directory if it doesn't exist
  my $debug_path = $self->conf->param('dumppath').'/debug';
  unless (-d $debug_path) {
    system("mkdir -p $debug_path") == 0 or
      throw("Unable to create directory $debug_path.\n");
  }
  
  my $logfile = "$debug_path/ambiguous_${type}.txt";
  
  open(my $fh, '>', $logfile) or
    throw("Unable to open $logfile for writing: $!");

  my @low_scoring = ();
  my @high_scoring = ();
  my $last_id;

  # log by source
  foreach my $entry (sort { $a->source <=> $b->source }
                        @{ $matrix->get_all_Entries }) {
    
    $last_id ||= $entry->target;

    if ($last_id != $entry->source) {
      $self->write_ambiguous($type, 'source', $fh, \@low_scoring,
        \@high_scoring);
      $last_id = $entry->source;
    }
    
    if ($entry->score < 0.5) {
      push @low_scoring, $entry;
    } else {
      push @high_scoring, $entry;
    }
  }

  # write last source
  $self->write_ambiguous($type, 'source', $fh, \@low_scoring, \@high_scoring);

  # now do the same by target
  $last_id = undef;
  foreach my $entry (sort { $a->target <=> $b->target }
                        @{ $matrix->get_all_Entries }) {

    $last_id ||= $entry->target;

    if ($last_id != $entry->target) {
      $self->write_ambiguous($type, 'target', $fh, \@low_scoring,
        \@high_scoring);
      $last_id = $entry->target;
    }
    
    if ($entry->score < 0.5) {
      push @low_scoring, $entry;
    } else {
      push @high_scoring, $entry;
    }
  }

  # write last target
  $self->write_ambiguous($type, 'target', $fh, \@low_scoring, \@high_scoring);

  close($fh);
}


sub write_ambiguous {
  my ($self, $type, $db_type, $fh, $low, $high) = @_;
  
  # if only source or target are ambiguous (i.e. you have only one mapping from
  # this perspective) then log from the other perspective
  if (scalar(@$low) + scalar(@$high) <= 1) {
    @$low = ();
    @$high = ();
    return;
  }

  my $first_id;
  if (@$low) {
    $first_id = $low->[0]->$db_type;
  } else {
    $first_id = $high->[0]->$db_type;
  }

  my $other_db_type;
  if ($db_type eq 'source') {
    $other_db_type = 'target';
  } else {
    $other_db_type = 'source';
  }

  print $fh "$db_type $type $first_id scores ambiguously:\n";

  # high scorers
  if (@$high) {
    print $fh "  high scoring ${other_db_type}s\n";

    while (my $e = shift(@$high)) {
      print $fh "    ", $e->$other_db_type, " ", $e->score, "\n";
    }
  }

  # low scorers
  if (@$low) {
    print $fh "  low scoring ${other_db_type}s\n    ";

    my $i = 1;

    while (my $e = shift(@$low)) {
      print $fh "\n    " unless (($i++)%10);
      print $fh $e->$other_db_type, ", ";
    }
  }

  print $fh "\n";
}


1;

