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

Bio::EnsEMBL::RepeatFeature - A feature representing a repeat on a piece of
sequence.

=head1 SYNOPSIS

  my $rf = new Bio::EnsEMBL::Feature(
    -start            => 100,
    -end              => 220,
    -strand           => -1,
    -slice            => $slice,
    -analysis         => $analysis,
    -repeat_consensus => $rc,
    -hstart           => 10,
    -hend             => 100,
    -hstrand          => 1,
    -score            => 83.2
  );

  my $hstart = $feat->hstart;
  my $hend   = $feat->hend;

  # move the feature to the chromosomal coordinate system
  $feature = $feature->transform('chromosome');

  # move the feature to a different slice
  # (possibly on another coord system)
  $feature = $feature->transfer($new_slice);

  # project the feature onto another coordinate system possibly across
  # boundaries:
  @projection = @{ $feature->project('contig') };

  # change the start, end, and strand of the feature in place
  $feature->move( $new_start, $new_end, $new_strand );

=head1 DESCRIPTION

This a feature representing a repeat region on a sequence

=head1 METHODS

=cut

package Bio::EnsEMBL::RepeatFeature;

use strict;

use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

use base qw/Bio::EnsEMBL::Feature/;

=head2 new

  Arg [REPEAT_CONSENSUS] : Bio::EnsEMBL::RepeatConsensus (optional)
                           The repeat consensus for this repeat feature
  Arg [HSTART] : int (optional)
                 The hit start on the consensus sequence
  Arg [HEND]   : int (optional)
                 The hit end on the consensus sequence
  Arg [SCORE]  : float (optional)
                 The score
  Arg [...]    : Named arguments to superclass constructor
                 (see Bio::EnsEMBL::Feaure)
  Example    : $rf = Bio::EnsEMBL::RepeatFeature->new(-REPEAT_CONSENSUS => $rc,
                                                      -HSTART => 10,
                                                      -HEND   => 100,
                                                      -SCORE  => 58.0,
                                                      -START  => 1_000_100,
                                                      -END    => 1_000_190,
                                                      -STRAND => 1,
                                                      -ANALYSIS => $an,
                                                      -SLICE  => $chr_slice);
  Description: Creates a new Bio::EnsEMBL::RepeatFeature object
  Returntype : Bio::EnsEMBL::RepeatFeature
  Exceptions : none
  Caller     : RepeatFeatureAdaptors
  Status     : Stable

=cut

sub new {
  my $caller = shift;

  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  my ($repeat_consensus, $hstart, $hend, $score) =
    rearrange(['REPEAT_CONSENSUS','HSTART','HEND','SCORE'], @_);

  $self->repeat_consensus($repeat_consensus);
  $self->{'hstart'} = $hstart;
  $self->{'hend'}   = $hend;
  $self->{'score'}  = $score;

  return $self;
}


=head2 repeat_consensus

  Arg [1]    : (optional) Bio::EnsEMBL::RepeatConsensus
  Example    : $repeat_consensus = $repeat->repeat_consensus;
  Description: Getter/Setter for the repeat consensus of this repeat
  Returntype : Bio::EnsEMBL::RepeatConsensus
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub repeat_consensus {
  my $self = shift;

  if(@_) {
    my $rc = shift;
    if(defined($rc)) {
      if(!ref($rc) || !$rc->isa('Bio::EnsEMBL::RepeatConsensus')) {
        throw('RepeatConsensus arg must be a Bio::EnsEMBL::RepeatConsensus');
      }
    }
    $self->{'repeat_consensus'} = $rc;
  }

  return $self->{'repeat_consensus'};
}



=head2 hstart

  Arg [1]    : (optional) int $hstart
  Example    : $hit_start = $repeat->hstart;
  Description: Getter/Setter for the start bp of this repeat match on the 
               consensus sequence.
  Returntype : int
  Exceptions : none 
  Caller     : general
  Status     : Stable

=cut

sub hstart {
  my $self = shift;
  $self->{'hstart'} = shift if(@_);
  return $self->{'hstart'};
}


=head2 score

  Arg [1]    : (optional) float $score
  Example    : $score = $repeat->score();
  Description: Getter/Setter for the score of this repeat feature
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub score {
  my $self = shift;
  $self->{'score'} = shift if(@_);
  return $self->{'score'};
}



=head2 hend

  Arg [1]    : (optional) int $hend
  Example    : $hit_end = $repeat->hend;
  Description: Getter/Setter for the end bp of this repeat match on the 
               consensus sequence.
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub hend {
  my $self = shift;
  $self->{'hend'} = shift if(@_);
  return $self->{'hend'};
}



=head2 hstrand

  Arg [1]    : none
  Example    : none
  Description: always returns 1. method exists for consistancy with other 
               features.
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub hstrand {
  return 1;
}


=head2 display_id

  Arg [1]    : none
  Example    : print $rf->display_id();
  Description: This method returns a string that is considered to be
               the 'display' identifier.  For repeat_features this is the 
               name of the repeat consensus if it is available otherwise it is
               an empty string.
  Returntype : string
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub display_id {
  my $self = shift;

  my $id = '';

  my $rc = $self->{'repeat_consensus'};
  if($rc) {
    $id = $rc->name();
  }

  return $id;
}


=head2 summary_as_hash

  Example       : $repeat_feature_summary = $protein_feature->summary_as_hash();
  Description   : Retrieves a textual summary of this Repeat feature.
                  Not inherited from Feature.
  Returns       : hashref of arrays of descriptive strings
  Status        : Intended for internal use
=cut

sub summary_as_hash {
  my $self = shift;
  my %summary;
  $summary{'start'} = $self->seq_region_start;
  $summary{'end'} = $self->seq_region_end;
  $summary{'strand'} = $self->strand;
  $summary{'seq_region_name'} = $self->seq_region_name;
  $summary{'description'} = $self->display_id;
  $summary{'assembly_name'} = $self->slice->coord_system->version();
  return \%summary;
}



1;

__END__

=head1 NAME - Bio::EnsEMBL::RepeatFeature

