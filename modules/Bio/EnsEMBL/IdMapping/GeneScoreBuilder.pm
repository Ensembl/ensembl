package Bio::EnsEMBL::IdMapping::GeneScoreBuilder;

=head1 NAME


=head1 SYNOPSIS


=head1 DESCRIPTION

Combines ExonScoreBuilder, ExonDirectMapper and ExonerateRunner from Java
application.

=head1 METHODS


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

use Bio::EnsEMBL::IdMapping::ScoreBuilder;
our @ISA = qw(Bio::EnsEMBL::IdMapping::ScoreBuilder);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::ScriptUtils qw(parse_bytes);
use Bio::EnsEMBL::IdMapping::ScoredMappingMatrix;


sub score_genes {
  my $self = shift;
  my $transcript_matrix = shift;

  unless ($transcript_matrix and
          $transcript_matrix->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix')) {
    throw('Need a Bio::EnsEMBL::IdMapping::ScoredMappingMatrix.');
  }

  $self->logger->log_stamped("Starting gene scoring...\n\n");

  # build scores based on transcript scores
  my $matrix = $self->scores_from_transcript_scores($transcript_matrix);

  # log stats of combined matrix
  my $fmt = "%-40s%10.0f\n";

  $self->logger->log("Scoring matrix:\n");

  $self->logger->log(sprintf($fmt, "Total source genes:",
    $self->cache->get_count_by_name('genes_by_id', 'source'), 1);

  $self->logger->log(sprintf($fmt, "Scored source genes:",
    $matrix->get_source_count), 1);

  $self->logger->log(sprintf($fmt, "Total target genes:",
    $self->cache->get_count_by_name('genes_by_id', 'target'), 1);

  $self->logger->log(sprintf($fmt, "Scored target genes:",
    $matrix->get_target_count), 1);

  $self->log_matrix_stats($matrix);
  
  $self->logger->log("\nDone with transcript scoring.\n\n");

  return $matrix;
}


sub scores_from_transcript_scores {
  my $self = shift;
  my $transcript_matrix = shift;

  unless ($transcript_matrix and
          $transcript_matrix->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix')) {
    throw('Need a Bio::EnsEMBL::IdMapping::ScoredMappingMatrix.');
  }
  
  my $matrix = Bio::EnsEMBL::IdMapping::ScoredMappingMatrix->new(
    -DUMP_PATH   => $self->conf->param('dumppath'),
    -CACHE_FILE  => 'gene_matrix.ser',
  );

  my $gene_cache = $matrix->cache_file;

  if (-s $gene_cache) {
    
    # read from file
    $self->logger->log_stamped("Reading gene scoring matrix from file...\n");
    $self->logger->log("Cache file $gene_cache.\n", 1);
    $matrix->read_from_file;
    $self->logger->log_stamped("Done.\n");
    
  } else {
    
    # build scoring matrix
    $self->logger->log("No gene scoring matrix found. Will build new one.\n");

    $self->logger->log_stamped("Transcript scoring...\n");
    $matrix = $self->build_scores($matrix, $transcript_matrix);
    $self->logger->log_stamped("Done.\n");

    # write scoring matrix to file
    $matrix->write_to_file;

  }

  return $matrix;
}


sub build_scores {
  my $self = shift;
  my $matrix = shift;
  my $transcript_matrix = shift;

  unless ($matrix and
          $matrix->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix')) {
    throw('Need a gene Bio::EnsEMBL::IdMapping::ScoredMappingMatrix.');
  }
  
  unless ($transcript_matrix and
          $transcript_matrix->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix')) {
    throw('Need a transcript Bio::EnsEMBL::IdMapping::ScoredMappingMatrix.');
  }

  # first find out which source and target genes have scoring transcripts and
  # build a "flag" matrix for these genes (all scores are 1)
  $self->flag_matrix_from_transcript_scores($matrix, $transcript_matrix);

  # now calculate the actual scores for the genes in the flag matrix
  $final_matrix = $self->score_matrix_from_flag_matrix($matrix,
    $transcript_matrix);
  
  return $final_matrix;
}


sub flag_matrix_from_transcript_scores {
  my $self = shift;
  my $matrix = shift;
  my $transcript_matrix = shift;

  # initialise progress logger
  my $i;
  my $num_genes =
    scalar(keys %{ $self->cache->get_by_name('genes_by_id', 'source') });
  
  $self->logger->log("Creating flag matrix...\n", 1);

  # for every transcript scoring matrix entry, make an entry in the gene flag
  # matrix.
  foreach my $entry (@{ $transcript_matrix->get_all_Entries }) {

    $self->logger->log_progress($num_genes, ++$i, 20, 1, 0);

    my $source_gene = $self->cache->get_by_key('genes_by_transcript_id',
      'source', $entry->source);

    my $target_gene = $self->cache->get_by_key('genes_by_transcript_id',
      'target', $entry->target);

    $matrix->add_score($source_gene->id, $target_gene->id, 1);
  }

  $self->logger->log("\n\n");

  return $matrix;
}


sub score_matrix_from_flag_matrix {
  my $self = shift;
  my $flag_matrix = shift;
  my $transcript_matrix = shift;
  my $simple_scoring = shift;

  # create a new scoring matrix which will replace the flag matrix
  my $matrix = Bio::EnsEMBL::IdMapping::ScoredMappingMatrix->new(
    -DUMP_PATH   => $self->conf->param('dumppath'),
    -CACHE_FILE  => 'gene_matrix.ser',
  );

  # initialise progress logger
  my $i;
  my $num_genes =
    scalar(keys %{ $self->cache->get_by_name('genes_by_id', 'source') });
  
  $self->logger->log("Creating score matrix from flag matrix...\n", 1);

  # loop over flag matrix and do proper scoring for each entry
  foreach my $entry (@{ $flag_matrix->get_all_Entries }) {
    
    $self->logger->log_progress($num_genes, ++$i, 20, 1, 0);

    my $score = 0;
    my $source_gene = $self->cache->get_by_key('genes_by_id', 'source',
      $entry->source);
    my $target_gene = $self->cache->get_by_key('genes_by_id', 'target',
      $entry->target);

    if ($simple_scoring) {
    
      # simple scoring (used for rescoring purposes
      $score = $self->simple_gene_gene_score($source_gene, $target_gene,
        $transcript_matrix);
    
    } else {

      # full scoring
      $score = $self->complex_gene_gene_score($source_gene, $target_gene,
        $transcript_matrix);
    }

    $matrix->add($entry->source, $entry->target, $score);
  }

  $self->logger->log("\n\n");

  return $matrix;
}


sub complex_gene_gene_score {
  my $self = shift;
  my $source_gene = shift;
  my $target_gene = shift;
  my $transcript_matrix = shift;

  my $source_gene_score = 0;
  my $target_gene_score = 0;
  my $source_gene_accum_length = 0; # sum of all transcript lengths
  my $target_gene_accum_length = 0; # sum of all transcript lengths

  # We are only interested in scoring with transcripts that are in the target
  # gene. The scored mapping matrix may contain scores for transcripts that
  # aren't in this gene so create a hash of the target genes's transcripts
  my %target_transcripts = map { $_->id => 1 }
    @{ $target_gene->get_all_Transcripts };
    
  # loop over source transcripts
  foreach my $source_transcript (@{ $source_gene->get_all_Transcripts ) {

    # now loop over target transcripts and find the highest scoring target
    # transcript belonging to the target gene
    my $max_source_score = -1;
    
    foreach my $target_transcript_id (@{ $transcript_matrix->get_targets_for_source($source_transcript->id) }) {

      next unless (%target_transcripts{$target_transcript_id});

      my $score = $transcript_matrix->get_score(
        $source_transcript->id, $target_transcript_id);
      $max_source_score = $score if ($score > $max_source_score);
    }

    if ($max_source_score > 0) {
      $source_gene_score += $max_source_score * $source_transcript->length;
    }

    $source_gene_accum_length += $source_transcript->length;
  }
    
  # now do the same for target genes
  my %source_transcripts = map { $_->id => 1 }
    @{ $source_gene->get_all_Transcripts };
    
  # loop over target transcripts
  foreach my $target_transcript (@{ $target_gene->get_all_Transcripts ) {

    # now loop over source transcripts and find the highest scoring source
    # transcript belonging to the source gene
    my $max_target_score = -1;
    
    foreach my $source_transcript_id (@{ $transcript_matrix->get_sources_for_target($target_transcript->id) }) {

      next unless (%source_transcripts{$source_transcript_id});

      my $score = $transcript_matrix->get_score(
        $source_transcript_id, $target_transcript->id);
      $max_target_score = $score if ($score > $max_target_score);
    }

    if ($max_target_score > 0) {
      $target_gene_score += $max_target_score * $target_transcript->length;
    }

    $target_gene_accum_length += $target_transcript->length;
  }

  # calculate overall score for this gene
  my $gene_score = 0;

  if (($source_gene_length + $target_gene_length) > 0) {

    $gene_score = ($source_gene_score + $target_gene_score) /
                  ($source_gene_length + $target_gene_length);

  } else {
  
    $self->logger->log_warning("Combined length of source (".$source_gene->id.") and target (".$target_gene->id.") gene is zero!\n", 1);
  
  }

  return $gene_score;
}


#
# Simplified scoring for genes. Score is best scoring transcript pair.
# This is used when the more elaborate gene representing score does not
# distinguish very well.
#
sub complex_gene_gene_score {
  my $self = shift;
  my $source_gene = shift;
  my $target_gene = shift;
  my $transcript_matrix = shift;

  my $gene_score = 0;

  foreach my $source_transcript (@{ $source_gene->get_all_Transcripts ) {
    foreach my $target_transcript (@{ $target_gene->get_all_Transcripts ) {
      
      my $score = $transcript_matrix->get_score($source_transcript->id,
        $target_transcript->id);

      $gene_score = $score if ($score > $gene_score);
    }
  }

  return $gene_score;
}


sub simple_gene_rescore {
  my $self = shift;
  my $gene_matrix = shift;
  my $transcript_matrix = shift;

  return $self->score_matrix_from_flag_matrix($gene_matrix,
    $transcript_matrix, 1);
}

1;

