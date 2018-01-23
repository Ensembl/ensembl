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

package Bio::EnsEMBL::IdMapping::TranscriptScoreBuilder;

use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::IdMapping::ScoreBuilder;
our @ISA = qw(Bio::EnsEMBL::IdMapping::ScoreBuilder);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::ScriptUtils qw(path_append);
use Bio::EnsEMBL::IdMapping::ScoredMappingMatrix;


sub score_transcripts {
  my $self = shift;
  my $exon_matrix = shift;

  unless ($exon_matrix and
          $exon_matrix->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix')) {
    throw('Need a Bio::EnsEMBL::IdMapping::ScoredMappingMatrix.');
  }

  $self->logger->info("-- Scoring transcripts...\n\n", 0, 'stamped');

  # build scores based on exon scores
  my $matrix = $self->scores_from_exon_scores($exon_matrix);

  # debug logging
  if ($self->logger->loglevel eq 'debug') {
    $matrix->log('transcript', $self->conf->param('basedir'));
  }

  # log stats of combined matrix
  my $fmt = "%-40s%10.0f\n";

  $self->logger->info("Scoring matrix:\n");

  $self->logger->info(sprintf($fmt, "Total source transcripts:",
    $self->cache->get_count_by_name('transcripts_by_id', 'source')), 1);

  $self->logger->info(sprintf($fmt, "Total target transcripts:",
    $self->cache->get_count_by_name('transcripts_by_id', 'target')), 1);

  $self->log_matrix_stats($matrix);
  
  $self->logger->info("\nDone with transcript scoring.\n\n");

  return $matrix;
}


sub scores_from_exon_scores {
  my $self = shift;
  my $exon_matrix = shift;

  unless ($exon_matrix and
          $exon_matrix->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix')) {
    throw('Need a Bio::EnsEMBL::IdMapping::ScoredMappingMatrix.');
  }
  
  my $dump_path = path_append($self->conf->param('basedir'), 'matrix');
  
  my $matrix = Bio::EnsEMBL::IdMapping::ScoredMappingMatrix->new(
    -DUMP_PATH   => $dump_path,
    -CACHE_FILE  => 'transcript_matrix.ser',
  );

  my $transcript_cache = $matrix->cache_file;

  if (-s $transcript_cache) {
    
    # read from file
    $self->logger->info("Reading transcript scoring matrix from file...\n", 0, 'stamped');
    $self->logger->debug("Cache file $transcript_cache.\n", 1);
    $matrix->read_from_file;
    $self->logger->info("Done.\n\n", 0, 'stamped');
    
  } else {
    
    # build scoring matrix
    $self->logger->info("No transcript scoring matrix found. Will build new one.\n");

    $self->logger->info("Transcript scoring...\n", 0, 'stamped');
    $matrix = $self->build_scores($matrix, $exon_matrix);
    $self->logger->info("Done.\n\n", 0, 'stamped');

    # write scoring matrix to file
    $matrix->write_to_file;

  }

  return $matrix;
}


sub build_scores {
  my $self = shift;
  my $matrix = shift;
  my $exon_matrix = shift;

  unless ($matrix and
          $matrix->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix')) {
    throw('Need a transcript Bio::EnsEMBL::IdMapping::ScoredMappingMatrix.');
  }
  
  unless ($exon_matrix and
          $exon_matrix->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix')) {
    throw('Need a exon Bio::EnsEMBL::IdMapping::ScoredMappingMatrix.');
  }

  # first find out which source and target transcripts have scoring exons and
  # build a "flag" matrix for these transcripts (all scores are 1)
  $self->flag_matrix_from_exon_scores($matrix, $exon_matrix);

  # now calculate the actual scores for the transcripts in the flag matrix
  my $final_matrix =
    $self->score_matrix_from_flag_matrix($matrix, $exon_matrix);
  
  return $final_matrix;
}


sub flag_matrix_from_exon_scores {
  my $self = shift;
  my $matrix = shift;
  my $exon_matrix = shift;

  # initialise progress logger
  my $i;
  my $num_transcripts =
    scalar(keys %{ $self->cache->get_by_name('transcripts_by_id', 'source') });
  my $progress_id = $self->logger->init_progress($num_transcripts, 100);

  $self->logger->info("Creating flag matrix...\n", 1);

  # loop over source transcripts
  foreach my $source_transcript (values %{ $self->cache->get_by_name('transcripts_by_id', 'source') }) {
    
    # log progress
    $self->logger->log_progress($progress_id, ++$i, 1);

    # get all exons for the source transcript
    foreach my $source_exon (@{ $source_transcript->get_all_Exons }) {

      # get target exons for this source exon from scoring matrix
      foreach my $target_exon_id (@{ $exon_matrix->get_targets_for_source($source_exon->id) }) {

        # get target transcripts that contain this exon
        foreach my $target_transcript (@{ $self->cache->get_by_key('transcripts_by_exon_id', 'target', $target_exon_id) }) {
          
          # add scoring flag for these two transcripts
          $matrix->add_score($source_transcript->id, $target_transcript->id, 1);
          
        }
      }
    }
  }

  $self->logger->info("\n");

  return $matrix;
}


sub score_matrix_from_flag_matrix {
  my $self = shift;
  my $flag_matrix = shift;
  my $exon_matrix = shift;

  unless ($flag_matrix and
          $flag_matrix->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix')) {
    throw('Need a transcript Bio::EnsEMBL::IdMapping::ScoredMappingMatrix.');
  }
  
  unless ($exon_matrix and
          $exon_matrix->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix')) {
    throw('Need an exon Bio::EnsEMBL::IdMapping::ScoredMappingMatrix.');
  }

  my $transcript_score_threshold =
    $self->conf->param('transcript_score_threshold') || 0;

  # create a new scoring matrix which will replace the flag matrix
  my $matrix = Bio::EnsEMBL::IdMapping::ScoredMappingMatrix->new(
    -DUMP_PATH   => $flag_matrix->dump_path,
    -CACHE_FILE  => $flag_matrix->cache_file_name,
  );

  # initialise progress logger
  my $i;
  my $num_transcripts =
    scalar(keys %{ $self->cache->get_by_name('transcripts_by_id', 'source') });
  my $progress_id = $self->logger->init_progress($num_transcripts, 100);

  $self->logger->info("Creating score matrix from flag matrix...\n", 1);

  # debug
  my $fmt_d1 = "%-14s%-15s%-14s%-14s%-14s\n";
  my $fmt_d2 = "%.6f";
  
  # loop over source transcripts
  foreach my $source_transcript (values %{ $self->cache->get_by_name('transcripts_by_id', 'source') }) {
    
    # log progress
    $self->logger->log_progress($progress_id, ++$i, 1);

    # We are only interested in scoring with exons that are in the target
    # transcript. The ScoredMappingMatrix may contain scores for exons that
    # aren't in this transcript so create a hash of the target transcript's
    # exons
    my %source_exons = map { $_->id => 1 }
      @{ $source_transcript->get_all_Exons };

    my $source_transcript_length = $source_transcript->length;

    # get all corresponding target transcripts from the flag matrix
    foreach my $target_transcript_id (@{ $flag_matrix->get_targets_for_source($source_transcript->id) }) {
      
      my $target_transcript = $self->cache->get_by_key('transcripts_by_id', 'target', $target_transcript_id);

      my $source_transcript_score = 0;
      my $target_transcript_score = 0;
      my $target_transcript_length = $target_transcript->length;

      my %target_exons = map { $_->id => 1 }
        @{ $target_transcript->get_all_Exons };

      # now loop over source exons and find the highest scoring target exon
      # belonging to the target transcript
      
      foreach my $source_exon (@{ $source_transcript->get_all_Exons }) {

        my $max_source_score = -1;
        
        foreach my $target_exon_id (@{ $exon_matrix->get_targets_for_source($source_exon->id) }) {

          next unless ($target_exons{$target_exon_id});

          my $score = $exon_matrix->get_score($source_exon->id,
            $target_exon_id);
          $max_source_score = $score if ($score > $max_source_score);
        }

        if ($max_source_score > 0) {
          $source_transcript_score += ($max_source_score * $source_exon->length);

        }
      }

      # now do the same for target exons
      
      foreach my $target_exon (@{ $target_transcript->get_all_Exons }) {

        my $max_target_score = -1;
        
        foreach my $source_exon_id (@{ $exon_matrix->get_sources_for_target($target_exon->id) }) {

          next unless ($source_exons{$source_exon_id});

          my $score = $exon_matrix->get_score(
            $source_exon_id, $target_exon->id);
          $max_target_score = $score if ($score > $max_target_score);
        }

        if ($max_target_score > 0) {
          $target_transcript_score += ($max_target_score * $target_exon->length);

        }
      }

      #
      # calculate transcript score and add to scoring matrix
      #
      if (($source_transcript_length + $target_transcript_length) > 0) {

        # sanity check
        if (($source_transcript_score > $source_transcript_length) or
            ($target_transcript_score > $target_transcript_length)) {

          $self->logger->warning("Score > length for source ($source_transcript_score <> $source_transcript_length) or target ($target_transcript_score <> $target_transcript_length).\n", 1);

        } else {

          my $source_transcript_biotype_group = $self->get_biotype_group($source_transcript->biotype());
          my $target_transcript_biotype_group = $self->get_biotype_group($target_transcript->biotype());

=cut
          # debug
          $self->logger->info($source_transcript->id.":".$target_transcript->id.
            " source score: $source_transcript_score".
            " source length: $source_transcript_length".
            " source biotype:" . $source_transcript->biotype() .
            " source group: $source_transcript_biotype_group".
            " target score: $target_transcript_score".
            " target biotype:" . $target_transcript->biotype() .
            " target group: $target_transcript_biotype_group".
            " target length: $target_transcript_length\n");
=cut
          
          my $transcript_score =
            ($source_transcript_score + $target_transcript_score) /
            ($source_transcript_length + $target_transcript_length);

## Add penalty if biotypes are different
          if ($source_transcript->biotype() ne $target_transcript->biotype()) {
            $transcript_score = $transcript_score * 0.9;
          }

## Add penalty if biotype groups are different
          if ($source_transcript_biotype_group ne $target_transcript_biotype_group) {
            $transcript_score = $transcript_score * 0.8;
          }

          # everything is fine, add score to matrix
          if ($transcript_score > $transcript_score_threshold) {
            $matrix->add_score($source_transcript->id, $target_transcript->id,
              $transcript_score);
          }
          
        }
      
      } else {
      
        $self->logger->warning("Combined length of source (".$source_transcript->id.") and target (".$target_transcript->id.") transcript is zero!\n", 1);
      
      }

    }
  }

  $self->logger->info("\n");

  return $matrix;
    
}


#
# penalise scores between genes with different biotypes.
# entries are modified in place
#
sub biotype_transcript_rescore {
  my ( $self, $matrix ) = @_;

  if ( defined($matrix) &&
       !$matrix->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix') )
  {
    throw('Exprected Bio::EnsEMBL::IdMapping::ScoredMappingMatrix');
  }

  my $i = 0;

  foreach my $entry ( @{ $matrix->get_all_Entries } ) {

    my $source_transcript =
      $self->cache->get_by_key( 'transcripts_by_id', 'source',
                                $entry->source() );

    my $target_transcript =
      $self->cache->get_by_key( 'transcripts_by_id', 'target',
                                $entry->target() );

    if ($source_transcript->biotype() ne $target_transcript->biotype() )
    {
      # PENALTY: Lower the score for mappings to transcripts of
      # different biotype.
      $matrix->set_score( $entry->source(), $entry->target(),
                          0.9*$entry->score() );
      $i++;
    }
  }

  $self->logger->debug("Scored transcripts with biotype mismatch: $i\n",
                       1 );
} ## end sub biotype_transcript_rescore


sub different_translation_rescore {
  my $self   = shift;
  my $matrix = shift;

  unless ($matrix
      and $matrix->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix') )
  {
    throw('Need a Bio::EnsEMBL::IdMapping::ScoredMappingMatrix.');
  }

  my $i = 0;

  foreach my $entry ( sort { $b->score <=> $a->score }
                      @{ $matrix->get_all_Entries } )
  {

    # we only do this for perfect matches, i.e. transcript score == 1
    last if ( $entry->score < 1 );

    my $source_tl =
      $self->cache->get_by_key( 'transcripts_by_id', 'source',
                                $entry->source )->translation;
    my $target_tl =
      $self->cache->get_by_key( 'transcripts_by_id', 'target',
                                $entry->target )->translation;

    # no penalty if both transcripts have no translation
    next if ( !$source_tl and !$target_tl );

    if (    !$source_tl
         or !$target_tl
         or ( $source_tl->seq ne $target_tl->seq ) )
    {
      # PENALTY: The transcript stable ID is now on a transcript with a
      # different translation.
      $matrix->set_score( $entry->source(), $entry->target(),
                          0.9*$entry->score() );
      $i++;
    }

  } ## end foreach my $entry ( sort { ...})

  $self->logger->debug(
                "Non-perfect translations on perfect transcripts: $i\n",
                1 );
} ## end sub different_translation_rescore


sub non_mapped_gene_rescore {
  my $self          = shift;
  my $matrix        = shift;
  my $gene_mappings = shift;

  # argument checks
  unless ($matrix
      and $matrix->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix') )
  {
    throw(
       'Need a transcript Bio::EnsEMBL::IdMapping::ScoredMappingMatrix.'
    );
  }

  unless ( $gene_mappings
       and $gene_mappings->isa('Bio::EnsEMBL::IdMapping::MappingList') )
  {
    throw('Need a gene Bio::EnsEMBL::IdMapping::MappingList.');
  }

  # create of lookup hash of mapped source genes to target genes
  my %gene_lookup =
    map { $_->source => $_->target }
    @{ $gene_mappings->get_all_Entries };

  my $i = 0;

  foreach my $entry ( @{ $matrix->get_all_Entries } ) {

    my $source_gene =
      $self->cache->get_by_key( 'genes_by_transcript_id', 'source',
                                $entry->source );
    my $target_gene =
      $self->cache->get_by_key( 'genes_by_transcript_id', 'target',
                                $entry->target );

    my $mapped_target = $gene_lookup{ $source_gene->id };

    if ( !$mapped_target or ( $mapped_target != $target_gene->id ) ) {
      # PENALTY: The transcript stable ID has been mapped to an
      # un-mapped gene.
      $matrix->set_score( $entry->source(), $entry->target(),
                          0.9*$entry->score() );
      $i++;
    }
  }

  $self->logger->debug( "Scored transcripts in non-mapped genes: $i\n",
                        1 );
} ## end sub non_mapped_gene_rescore


sub get_biotype_group {
  my ($self, $biotype) = @_;
  my $dba = $self->cache->get_production_DBAdaptor();
  my $helper   = $self->cache->get_production_DBAdaptor()->dbc()->sql_helper();

  my $sql      = q{
     SELECT biotype_group
     FROM biotype
     WHERE object_type = 'transcript'
     AND is_current = 1
     AND name = ?
     AND db_type like '%core%' };
  my $result = $helper->execute_simple(-SQL => $sql, -PARAMS => [$biotype]);
  return $result->[0];
}

  
1;

