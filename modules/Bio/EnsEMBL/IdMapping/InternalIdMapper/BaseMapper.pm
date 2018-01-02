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

package Bio::EnsEMBL::IdMapping::InternalIdMapper::BaseMapper;

use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::IdMapping::BaseObject;
our @ISA = qw(Bio::EnsEMBL::IdMapping::BaseObject);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::ScriptUtils qw(path_append);
use Bio::EnsEMBL::IdMapping::MappingList;

# scores are considered the same if (2.0 * (s1-s2))/(s1 + s2) < this
use constant SIMILAR_SCORE_RATIO => 0.01;

#
# find the highest unambiguous score for all sources and targets in a scoring
# matrix
#
sub basic_mapping {
  my $self         = shift;
  my $matrix       = shift;
  my $mapping_name = shift;

  # argument checks
  unless ($matrix
      and $matrix->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix') )
  {
    throw('Need a Bio::EnsEMBL::IdMapping::ScoredMappingMatrix.');
  }

  throw('Need a name for serialising the mapping.')
    unless ($mapping_name);

  # Create a new MappingList object. Specify AUTO_LOAD to load
  # serialised existing mappings if found
  my $dump_path =
    path_append( $self->conf->param('basedir'), 'mapping' );

  my $mappings =
    Bio::EnsEMBL::IdMapping::MappingList->new(
                                   -DUMP_PATH  => $dump_path,
                                   -CACHE_FILE => "${mapping_name}.ser",
                                   -AUTO_LOAD  => 1, );

  # checkpoint test: return a previously stored MappingList
  if ( $mappings->loaded ) {
    $self->logger->info(
                  "Read existing mappings from ${mapping_name}.ser.\n");
    return $mappings;
  }

  my $sources_done = {};
  my $targets_done = {};

  # sort scoring matrix entries by descending score
  my @sorted_entries =
    sort { $b->score <=> $a->score } @{ $matrix->get_all_Entries };

  # debug
  #my $idx = substr($mapping_name, -1);

  while ( my $entry = shift(@sorted_entries) ) {

    #$self->logger->debug("\nxxx$idx ".$entry->to_string." ");

    # we already found a mapping for either source or target
    next
      if (    $sources_done->{ $entry->source }
           or $targets_done->{ $entry->target } );

    #$self->logger->debug('d');

    # there's a better mapping for either source or target
    next
      if ( $self->higher_score_exists(
                           $entry, $matrix, $sources_done, $targets_done
           ) );

    #$self->logger->debug('h');

    # check for ambiguous mappings; they are dealt with later
    my $other_sources = [];
    my $other_targets = [];

    if ( $self->ambiguous_mapping( $entry,         $matrix,
                                   $other_sources, $other_targets ) )
    {
      #$self->logger->debug('a');

      $other_sources =
        $self->filter_sources( $other_sources, $sources_done );
      $other_targets =
        $self->filter_targets( $other_targets, $targets_done );

      next if ( scalar(@$other_sources) or scalar(@$other_targets) );
    }

    #$self->logger->debug('A');

    # this is the best mapping, add it
    $mappings->add_Entry($entry);

    $sources_done->{ $entry->source } = 1;
    $targets_done->{ $entry->target } = 1;
  } ## end while ( my $entry = shift...)

  # create checkpoint
  $mappings->write_to_file;

  return $mappings;
} ## end sub basic_mapping

sub higher_score_exists {
  my ( $self, $entry, $matrix, $sources_done, $targets_done ) = @_;

  my $source = $entry->source;
  my $target = $entry->target;
  my $score  = $entry->score;

  foreach
    my $other_source ( @{ $matrix->get_sources_for_target($target) } )
  {
    if (     $other_source != $source
         and !$sources_done->{$other_source}
         and $score < $matrix->get_score( $other_source, $target ) )
    {
      return 1;
    }
  }

  foreach
    my $other_target ( @{ $matrix->get_targets_for_source($source) } )
  {
    if (     $other_target != $target
         and !$targets_done->{$other_target}
         and $score < $matrix->get_score( $source, $other_target ) )
    {
      return 1;
    }
  }

  return 0;
} ## end sub higher_score_exists

#
# find ambiguous mappings (see scores_similar() for definition)
#
sub ambiguous_mapping {
  my ( $self, $entry, $matrix, $other_sources, $other_targets ) = @_;

  my $source = $entry->source;
  my $target = $entry->target;
  my $score  = $entry->score;

  my $retval = 0;

  foreach
    my $other_source ( @{ $matrix->get_sources_for_target($target) } )
  {
    my $other_score = $matrix->get_score( $other_source, $target );

    if ( $other_source != $source
         and (    $self->scores_similar( $score, $other_score )
               or $score < $other_score ) )
    {
      $retval = 1;
      push @{$other_sources}, $other_source;
    }
  }

  foreach
    my $other_target ( @{ $matrix->get_targets_for_source($source) } )
  {
    my $other_score = $matrix->get_score( $source, $other_target );

    if ( $other_target != $target
         and (    $self->scores_similar( $score, $other_score )
               or $score < $other_score ) )
    {
      $retval = 1;
      push @{$other_targets}, $other_target;
    }
  }

  return $retval;
} ## end sub ambiguous_mapping

#
# rule for similarity taken from java code...
#
sub scores_similar {
  my ( $self, $s1, $s2 ) = @_;

  # always give priority to exact matches over very similar ones
  return 0 if ( $s1 == 1 and $s2 < 1 );

  my $diff = $s1 - $s2;
  $diff = -$diff if ( $diff < 0 );

  my $pc = 2*$diff/( $s1 + $s2 );

  return ( $pc < SIMILAR_SCORE_RATIO );
}

sub filter_sources {
  my ( $self, $other_sources, $sources_done ) = @_;

  unless (     scalar( @{$other_sources} )
           and scalar( keys %{$sources_done} ) )
  {
    return $other_sources;
  }

  my @tmp = ();

  foreach my $e ( @{$other_sources} ) {
    push @tmp, $e unless ( $sources_done->{$e} );
  }

  return \@tmp;
}

sub filter_targets {
  my ( $self, $other_targets, $targets_done ) = @_;

  unless (     scalar( @{$other_targets} )
           and scalar( keys %{$targets_done} ) )
  {
    return $other_targets;
  }

  my @tmp = ();

  foreach my $e ( @{$other_targets} ) {
    push @tmp, $e unless ( $targets_done->{$e} );
  }

  return \@tmp;
}

1;
