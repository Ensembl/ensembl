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

=head1 METHODS

=cut


package Bio::EnsEMBL::IdMapping::InternalIdMapper::EnsemblTranscriptGeneric;

use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::IdMapping::InternalIdMapper::BaseMapper;
our @ISA = qw(Bio::EnsEMBL::IdMapping::InternalIdMapper::BaseMapper);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::ScriptUtils qw(path_append);

  
#
# basic mapping
#
sub init_basic {
  my $self = shift;
  my $num = shift;
  my $tsb = shift;
  my $mappings = shift;
  my $transcript_scores = shift;

  $self->logger->info("Basic transcript mapping...\n", 0, 'stamped');

  $mappings = $self->basic_mapping($transcript_scores,
    "transcript_mappings$num");
  $num++;
  my $new_scores = $tsb->create_shrinked_matrix($transcript_scores, $mappings,
    "transcript_matrix$num");

  return ($new_scores, $mappings);
}


#
# handle cases with exact match but different translation
#
sub non_exact_translation {
  my $self = shift;
  my $num = shift;
  my $tsb = shift;
  my $mappings = shift;
  my $transcript_scores = shift;

  $self->logger->info("Exact Transcript non-exact Translation...\n", 0, 'stamped');
  
  unless ($transcript_scores->loaded) {
    $tsb->different_translation_rescore($transcript_scores);
    $transcript_scores->write_to_file;
  }
  
  $mappings = $self->basic_mapping($transcript_scores,
    "transcript_mappings$num");
  $num++;
  my $new_scores = $tsb->create_shrinked_matrix($transcript_scores, $mappings,
    "transcript_matrix$num");

  return ($new_scores, $mappings);
}


#
# reduce score for mappings of transcripts which do not belong to mapped
# genes
#
sub mapped_gene {
  my $self = shift;
  my $num = shift;
  my $tsb = shift;
  my $mappings = shift;
  my $transcript_scores = shift;
  my $gene_mappings = shift;

  $self->logger->info("Transcripts in mapped genes...\n", 0, 'stamped');
  
  unless ($transcript_scores->loaded) {
  $tsb->non_mapped_gene_rescore($transcript_scores, $gene_mappings);
    $transcript_scores->write_to_file;
  }
  
  $mappings = $self->basic_mapping($transcript_scores,
    "transcript_mappings$num");
  $num++;
  my $new_scores = $tsb->create_shrinked_matrix($transcript_scores, $mappings,
    "transcript_matrix$num");

  return ($new_scores, $mappings);
}

#
# rescore by penalising scores between transcripts with different biotypes
#
sub biotype {
  my $self              = shift;
  my $num               = shift;
  my $tsb               = shift;
  my $mappings          = shift;
  my $transcript_scores = shift;

  $self->logger->info( "Retry with biotype disambiguation...\n",
                       0, 'stamped' );

  unless ( $transcript_scores->loaded() ) {
    $tsb->biotype_transcript_rescore($transcript_scores);
    $transcript_scores->write_to_file();
  }

  my $new_mappings = $self->basic_mapping( $transcript_scores,
                                           "transcript_mappings$num" );
  $num++;
  my $new_scores =
    $tsb->create_shrinked_matrix( $transcript_scores, $new_mappings,
                                  "transcript_matrix$num" );

  return ( $new_scores, $new_mappings );
}

#
# selectively rescore by penalising scores between transcripts with
# different internalIDs  
#
sub internal_id {
  my $self = shift;
  my $num = shift;
  my $tsb = shift;
  my $mappings = shift;
  my $transcript_scores = shift;

  $self->logger->info("Retry with internalID disambiguation...\n", 0, 'stamped');
  
  unless ($transcript_scores->loaded) {
    $tsb->internal_id_rescore($transcript_scores);
    $transcript_scores->write_to_file;
  }

  $mappings = $self->basic_mapping($transcript_scores,
    "transcript_mappings$num");
  $num++;
  my $new_scores = $tsb->create_shrinked_matrix($transcript_scores, $mappings,
    "transcript_matrix$num");

  return ($new_scores, $mappings);
}


#
# handle ambiguities between transcripts in single genes
#
sub single_gene {
  my $self = shift;
  my $num = shift;
  my $tsb = shift;
  my $mappings = shift;
  my $transcript_scores = shift;

  $self->logger->info("Transcripts in single genes...\n", 0, 'stamped');
  
  unless ($transcript_scores->loaded) {
    $transcript_scores->write_to_file;
  }
  
  $mappings = $self->same_gene_transcript_mapping($transcript_scores,
    "transcript_mappings$num");
  $num++;
  my $new_scores = $tsb->create_shrinked_matrix($transcript_scores, $mappings,
    "transcript_matrix$num");

  return ($new_scores, $mappings);
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
  my $dump_path = path_append($self->conf->param('basedir'), 'mapping');
  
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


1;

