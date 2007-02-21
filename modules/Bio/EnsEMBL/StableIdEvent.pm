package Bio::EnsEMBL::StableIdEvent;

=head1 NAME

Bio::EnsEMBL::StableIdEvent- object representing a stable ID mapping event

=head1 SYNOPSIS


=head1 DESCRIPTION

This object represents a stable ID mapping event. Such an event links two
ArchiveStableIds with a mapping score.

=head1 METHODS

score
old_ArchiveStableId
new_ArchiveStableId

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

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);


=head2 new

  Arg[1]      : Bio::EnsEMBL::ArchiveStableId $old_id
                The old ArchiveStableId in the mapping event
  Arg[2]      : Bio::EnsEMBL::ArchiveStableId $new_id
                The new ArchiveStableId in the mapping event
  Arg[3]      : (optional) float $score - score of this mapping event
  Example     : my $event = Bio::EnsEMBL::StableIdEvent->new(
                  $arch_id1, $arch_id2, 0.977);
  Description : object constructor
  Return type : Bio::EnsEMBL::StableIdEvent
  Exceptions  : thrown on wrong argument types
  Caller      : Bio::EnsEMBL::DBSQL::ArchiveStableIdAdaptor::fetch_history_tree_by_stable_id, general
  Status      : At Risk
              : under development

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($old_id, $new_id, $score) = rearrange([qw(OLD_ID NEW_ID SCORE)], @_);

  throw("Need old or new Bio::EnsEMBL::ArchiveStableId to create StableIdEvent")
    unless ($old_id || $new_id);

  my $self = {};
  bless $self, $class;

  # initialise object
  $self->old_ArchiveStableId($old_id);
  $self->new_ArchiveStableId($new_id);
  $self->score($score);
  
  return $self;
}


=head2 old_ArchiveStableId

  Arg[1]      : (optional) Bio::EnsEMBL::ArchiveStableId $archive_id, or undef
                The old ArchiveStableId to set for this mapping event
  Example     : # getter
                my $archive_id = $event->old_ArchiveStableId;
                
                # setter
                $event->old_ArchiveStableId($archive_id);
  Description : Getter/setter for old ArchiveStableId in this mapping event.
  Return type : Bio::EnsEMBL::ArchiveStableId
  Exceptions  : thrown on wrong argument type
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub old_ArchiveStableId {
  my $self = shift;
  
  # setter
  if (@_) {
    my $archive_id = shift;

    # if argument is defined, check type. undef is also legal as an argument.
    if (defined($archive_id)) {
      throw("Need a Bio::EnsEMBL::ArchiveStableId.") unless
        (ref($archive_id) && $archive_id->isa('Bio::EnsEMBL::ArchiveStableId'));
    }

    $self->{'old_id'} = $archive_id;
  }

  # getter
  return $self->{'old_id'};
}


=head2 new_ArchiveStableId

  Arg[1]      : (optional) Bio::EnsEMBL::ArchiveStableId $archive_id, or undef
                The new ArchiveStableId to set for this mapping event
  Example     : # getter
                my $archive_id = $event->new_ArchiveStableId;
                
                # setter
                $event->new_ArchiveStableId($archive_id);
  Description : Getter/setter for new ArchiveStableId in this mapping event.
  Return type : Bio::EnsEMBL::ArchiveStableId
  Exceptions  : thrown on wrong argument type
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub new_ArchiveStableId {
  my $self = shift;
  
  # setter
  if (@_) {
    my $archive_id = shift;

    # if argument is defined, check type. undef is also legal as an argument.
    if (defined($archive_id)) {
      throw("Need a Bio::EnsEMBL::ArchiveStableId.") unless
        (ref($archive_id) && $archive_id->isa('Bio::EnsEMBL::ArchiveStableId'));
    }

    $self->{'new_id'} = $archive_id;
  }

  # getter
  return $self->{'new_id'};
}


=head2 score

  Arg[1]      : (optional) float $score - the score to set
  Example     : my $score = $event->score;
  Description : Getter/setter for mapping event score.
  Return type : float or undef
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub score {
  my $self = shift;
  $self->{'score'} = shift if (@_);
  return $self->{'score'};
}


1;

