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

package Bio::EnsEMBL::IdMapping::InternalIdMapper::EnsemblExonGeneric;

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
  my $esb = shift;
  my $mappings = shift;
  my $exon_scores = shift;

  $self->logger->info("Basic exon mapping...\n", 0, 'stamped');

  $mappings = $self->basic_mapping($exon_scores, "exon_mappings$num");
  $num++;
  my $new_scores = $esb->create_shrinked_matrix($exon_scores, $mappings,
    "exon_matrix$num");

  return ($new_scores, $mappings);
}


#
# reduce score for mappings of exons which do not belong to mapped
# transcripts
# (ie where source exon transcript does not map target exon transcript)
#
sub mapped_transcript {
  my $self = shift;
  my $num = shift;
  my $esb = shift;
  my $mappings = shift;
  my $exon_scores = shift;
  my $transcript_mappings = shift;

  $self->logger->info("Exons in mapped transcript...\n", 0, 'stamped');

  unless ($exon_scores->loaded) {
    $esb->non_mapped_transcript_rescore($exon_scores, $transcript_mappings);
    $exon_scores->write_to_file;
  }

  $mappings = $self->basic_mapping($exon_scores, "exon_mappings$num");
  $num++;
  my $new_scores = $esb->create_shrinked_matrix($exon_scores, $mappings,
    "exon_matrix$num");

  return ($new_scores, $mappings);
}

sub single_transcript {
  my $self = shift;
  my $num = shift;
  my $esb = shift;
  my $mappings = shift;
  my $exon_scores = shift;

  $self->logger->info("Exons in single transcript...\n", 0, 'stamped');

  unless ($exon_scores->loaded) {
    $exon_scores->write_to_file;
  }

  $mappings = $self->same_transcript_exon_mapping($exon_scores, "exon_mappings$num");
  $num++;
  my $new_scores = $esb->create_shrinked_matrix($exon_scores, $mappings,
    "exon_matrix$num");

  return ($new_scores, $mappings);
}

sub same_transcript_exon_mapping {
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
    my %source_transcripts = ();
    my %target_transcripts = ();

    if ($self->ambiguous_mapping($entry, $matrix, $other_sources, $other_targets)) {
      #$self->logger->debug('a');

      $other_sources = $self->filter_sources($other_sources, $sources_done);
      $other_targets = $self->filter_targets($other_targets, $targets_done);

      $source_transcripts{$self->cache->get_by_key('transcripts_by_exon_id',
        'source', $entry->source)} = 1;
      $target_transcripts{$self->cache->get_by_key('transcripts_by_exon_id',
        'target', $entry->target)} = 1;

      foreach my $other_source (@{ $other_sources }) {
        $source_transcripts{$self->cache->get_by_key('transcripts_by_exon_id',
          'source', $other_source)} = 1;
      }

      foreach my $other_target (@{ $other_targets }) {
        $target_transcripts{$self->cache->get_by_key('transcripts_by_exon_id',
          'target', $other_target)} = 1;
      }

      # only add mapping if only one source and target gene involved
      if (scalar(keys %source_transcripts) == 1 and scalar(keys %target_transcripts) == 1) {
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
  

#
# selectively rescore by penalising scores between exons with
# different internalIDs
#
sub internal_id {
  my $self        = shift;
  my $num         = shift;
  my $esb         = shift;
  my $mappings    = shift;
  my $exon_scores = shift;

  $self->logger->info( "Retry with internalID disambiguation...\n",
                       0, 'stamped' );

  if ( !$exon_scores->loaded() ) {
    $esb->internal_id_rescore($exon_scores);
    $exon_scores->write_to_file();
  }

  $mappings = $self->basic_mapping( $exon_scores, "exon_mappings$num" );
  $num++;
  my $new_scores =
    $esb->create_shrinked_matrix( $exon_scores, $mappings,
                                  "exon_matrix$num" );

  return ( $new_scores, $mappings );
}


1;

