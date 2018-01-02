=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

=head1 SYNOPSIS

=head1 DESCRIPTION

Combines ExonScoreBuilder, ExonDirectMapper and ExonerateRunner from
Java application.

=head1 METHODS

=cut

package Bio::EnsEMBL::IdMapping::GeneScoreBuilder;

use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::IdMapping::ScoreBuilder;
our @ISA = qw(Bio::EnsEMBL::IdMapping::ScoreBuilder);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::ScriptUtils qw(path_append);
use Bio::EnsEMBL::IdMapping::ScoredMappingMatrix;


sub score_genes {
  my $self = shift;
  my $transcript_matrix = shift;

  unless ($transcript_matrix and
          $transcript_matrix->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix')) {
    throw('Need a Bio::EnsEMBL::IdMapping::ScoredMappingMatrix.');
  }

  $self->logger->info("-- Scoring genes...\n\n", 0, 'stamped');

  # build scores based on transcript scores
  my $matrix = $self->scores_from_transcript_scores($transcript_matrix);

  # debug logging
  if ($self->logger->loglevel eq 'debug') {
    $matrix->log('gene', $self->conf->param('basedir'));
  }

  # log stats of combined matrix
  my $fmt = "%-40s%10.0f\n";

  $self->logger->info("Scoring matrix:\n");

  $self->logger->info(sprintf($fmt, "Total source genes:",
    $self->cache->get_count_by_name('genes_by_id', 'source')), 1);

  $self->logger->info(sprintf($fmt, "Total target genes:",
    $self->cache->get_count_by_name('genes_by_id', 'target')), 1);

  $self->log_matrix_stats($matrix);
  
  $self->logger->info("\nDone with gene scoring.\n\n");

  return $matrix;
}


sub scores_from_transcript_scores {
  my $self = shift;
  my $transcript_matrix = shift;

  unless ($transcript_matrix and
          $transcript_matrix->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix')) {
    throw('Need a Bio::EnsEMBL::IdMapping::ScoredMappingMatrix.');
  }
  
  my $dump_path = path_append($self->conf->param('basedir'), 'matrix');
  
  my $matrix = Bio::EnsEMBL::IdMapping::ScoredMappingMatrix->new(
    -DUMP_PATH   => $dump_path,
    -CACHE_FILE  => 'gene_matrix.ser',
  );

  my $gene_cache = $matrix->cache_file;

  if (-s $gene_cache) {
    
    # read from file
    $self->logger->info("Reading gene scoring matrix from file...\n", 0, 'stamped');
    $self->logger->debug("Cache file $gene_cache.\n", 1);
    $matrix->read_from_file;
    $self->logger->info("Done.\n\n", 0, 'stamped');
    
  } else {
    
    # build scoring matrix
    $self->logger->info("No gene scoring matrix found. Will build new one.\n");

    $self->logger->info("Gene scoring...\n", 0, 'stamped');
    $matrix = $self->build_scores($matrix, $transcript_matrix);
    $self->logger->info("Done.\n\n", 0, 'stamped');

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
  my $final_matrix = $self->score_matrix_from_flag_matrix($matrix,
    $transcript_matrix);
  
  return $final_matrix;
}


sub flag_matrix_from_transcript_scores {
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

  # initialise progress logger
  my $i;
  my $num_entries = $transcript_matrix->get_entry_count;
  my $progress_id = $self->logger->init_progress($num_entries, 100);
  
  $self->logger->info("Creating flag matrix...\n", 1);

  # for every transcript scoring matrix entry, make an entry in the gene flag
  # matrix.
  foreach my $entry (@{ $transcript_matrix->get_all_Entries }) {

    $self->logger->log_progress($progress_id, ++$i, 1);

    my $source_gene = $self->cache->get_by_key('genes_by_transcript_id',
      'source', $entry->source);

    my $target_gene = $self->cache->get_by_key('genes_by_transcript_id',
      'target', $entry->target);

    $matrix->add_score($source_gene->id, $target_gene->id, 1);
  }

  $self->logger->info("\n");

  return $matrix;
}


sub score_matrix_from_flag_matrix {
  my $self              = shift;
  my $flag_matrix       = shift;
  my $transcript_matrix = shift;
  my $simple_scoring    = shift;

  unless ( $flag_matrix and
     $flag_matrix->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix') )
  {
    throw('Need a gene Bio::EnsEMBL::IdMapping::ScoredMappingMatrix.');
  }

  unless (     $transcript_matrix
           and
           $transcript_matrix->isa(
                         'Bio::EnsEMBL::IdMapping::ScoredMappingMatrix')
    )
  {
    throw(
       'Need a transcript Bio::EnsEMBL::IdMapping::ScoredMappingMatrix.'
    );
  }

  # create a new scoring matrix which will replace the flag matrix
  my $matrix =
    Bio::EnsEMBL::IdMapping::ScoredMappingMatrix->new(
                           -DUMP_PATH  => $flag_matrix->dump_path,
                           -CACHE_FILE => $flag_matrix->cache_file_name,
    );

  # initialise progress logger
  my $i;
  my $num_entries = $flag_matrix->get_entry_count;
  my $progress_id = $self->logger->init_progress( $num_entries, 100 );

  $self->logger->info( "Creating score matrix from flag matrix...\n",
                       1 );

  my $gene_score_threshold =
    $self->conf()->param('gene_score_threshold') || 0;

  # loop over flag matrix and do proper scoring for each entry
  foreach my $entry ( @{ $flag_matrix->get_all_Entries } ) {

    $self->logger->log_progress( $progress_id, ++$i, 1 );

    my $score = 0;
    my $source_gene =
      $self->cache->get_by_key( 'genes_by_id', 'source',
                                $entry->source );
    my $target_gene =
      $self->cache->get_by_key( 'genes_by_id', 'target',
                                $entry->target );

    if ($simple_scoring) {

      # simple scoring (used for rescoring purposes)
      $score =
        $self->simple_gene_gene_score( $source_gene, $target_gene,
                                       $transcript_matrix );

    }
    else {

      # full scoring
      $score =
        $self->complex_gene_gene_score( $source_gene, $target_gene,
                                        $transcript_matrix );
    }

    if ( $score > $gene_score_threshold ) {
      $matrix->add_score( $entry->source, $entry->target, $score );
    }
  } ## end foreach my $entry ( @{ $flag_matrix...})

  $self->logger->info("\n");

  return $matrix;
} ## end sub score_matrix_from_flag_matrix


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
  foreach my $source_transcript (@{ $source_gene->get_all_Transcripts }) {

    # now loop over target transcripts and find the highest scoring target
    # transcript belonging to the target gene
    my $max_source_score = -1;
    
    foreach my $target_transcript_id (@{ $transcript_matrix->get_targets_for_source($source_transcript->id) }) {

      next unless ($target_transcripts{$target_transcript_id});

      my $score = $transcript_matrix->get_score(
        $source_transcript->id, $target_transcript_id);
      $max_source_score = $score if ($score > $max_source_score);
    }

    if ($max_source_score > 0) {
      $source_gene_score += ($max_source_score * $source_transcript->length);
    }

    $source_gene_accum_length += $source_transcript->length;
  }
    
  # now do the same for target genes
  my %source_transcripts = map { $_->id => 1 }
    @{ $source_gene->get_all_Transcripts };
    
  # loop over target transcripts
  foreach my $target_transcript (@{ $target_gene->get_all_Transcripts }) {

    # now loop over source transcripts and find the highest scoring source
    # transcript belonging to the source gene
    my $max_target_score = -1;
    
    foreach my $source_transcript_id (@{ $transcript_matrix->get_sources_for_target($target_transcript->id) }) {

      next unless ($source_transcripts{$source_transcript_id});

      my $score = $transcript_matrix->get_score(
        $source_transcript_id, $target_transcript->id);
      $max_target_score = $score if ($score > $max_target_score);
    }

    if ($max_target_score > 0) {
      $target_gene_score += ($max_target_score * $target_transcript->length);
    }

    $target_gene_accum_length += $target_transcript->length;
  }

  # calculate overall score for this gene
  my $gene_score = 0;

  if (($source_gene_accum_length + $target_gene_accum_length) > 0) {

    $gene_score = ($source_gene_score + $target_gene_score) /
                  ($source_gene_accum_length + $target_gene_accum_length);

  } else {
  
    $self->logger->warning("Combined transcript length of source (".$source_gene->id.") and target (".$target_gene->id.") gene is zero!\n", 1);
  
  }

  if ($gene_score > 1) {
    $self->logger->warning("Illegal gene score: $gene_score (".
      join("|", $source_gene_score, $target_gene_score,
           $source_gene_accum_length, $target_gene_accum_length).")\n", 1);
  }

  return $gene_score;
}


#
# Simplified scoring for genes. Score is best scoring transcript pair.
# This is used when the more elaborate gene representing score does not
# distinguish very well.
#
sub simple_gene_gene_score {
  my $self = shift;
  my $source_gene = shift;
  my $target_gene = shift;
  my $transcript_matrix = shift;

  my $gene_score = 0;

  foreach my $source_transcript (@{ $source_gene->get_all_Transcripts }) {
    foreach my $target_transcript (@{ $target_gene->get_all_Transcripts }) {
      
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

  $gene_matrix = $self->score_matrix_from_flag_matrix($gene_matrix,
    $transcript_matrix, 1);
}

#
# penalise scores between genes with different biotypes.
# entries are modified in place
#
sub biotype_gene_rescore {
  my $self   = shift;
  my $matrix = shift;

  unless ($matrix
      and $matrix->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix') )
  {
    throw('Need a Bio::EnsEMBL::IdMapping::ScoredMappingMatrix.');
  }

  my $i = 0;

  foreach my $entry ( @{ $matrix->get_all_Entries } ) {

    my $source_gene =
      $self->cache->get_by_key( 'genes_by_id', 'source',
                                $entry->source );

    my $target_gene =
      $self->cache->get_by_key( 'genes_by_id', 'target',
                                $entry->target );

    if ( $source_gene->biotype() ne $target_gene->biotype() ) {
      # PENALTY: Lower the score for mappings that differ in biotype.
      $matrix->set_score( $entry->source(), $entry->target(),
                          0.9*$entry->score() );
      $i++;
    }
  }

  $self->logger->debug( "Scored genes with biotype mismatch: $i\n", 1 );
} ## end sub biotype_gene_rescore


sub location_gene_rescore {
  my $self   = shift;
  my $matrix = shift;

  unless ($matrix
      and $matrix->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix') )
  {
    throw('Need a Bio::EnsEMBL::IdMapping::ScoredMappingMatrix.');
  }

  my $i = 0;

  foreach my $entry ( @{ $matrix->get_all_Entries } ) {

    my $source_gene =
      $self->cache->get_by_key( 'genes_by_id', 'source',
                                $entry->source );

    my $target_gene =
      $self->cache->get_by_key( 'genes_by_id', 'target',
                                $entry->target );

    if ( $source_gene->seq_region_name() ne $target_gene->seq_region_name() ) {
      # PENALTY: Lower the score for mappings that are not on the same slice.
      $matrix->set_score( $entry->source(), $entry->target(),
                          0.9*$entry->score() );
      $i++;
    }
  }

}


1;

