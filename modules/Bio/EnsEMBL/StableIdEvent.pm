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

Bio::EnsEMBL::StableIdEvent- object representing a stable ID mapping event

=head1 SYNOPSIS

  my $old_id = Bio::EnsEMBL::ArchiveStableId->new(
    -stable_id => 'ENSG001',
    -version   => 1,
    -type      => 'Gene',
  );

  my $new_id = Bio::EnsEMBL::ArchiveStableId->new(
    -stable_id => 'ENSG001',
    -version   => 2,
    -type      => 'Gene',
  );

  my $event = Bio::EnsEMBL::StableIdEvent->new(
    -old_id => $old_id,
    -new_id => $new_id,
    -score  => 0.997
  );

  # directly access attributes in old and new ArchiveStableId
  my $old_stable_id = $event->get_attribute( 'old', 'stable_id' );

=head1 DESCRIPTION

This object represents a stable ID mapping event. Such an event links two
ArchiveStableIds with a mapping score.

=head1 METHODS

  new
  old_ArchiveStableId
  new_ArchiveStableId
  score
  get_attribute
  ident_string

=head1 RELATED MODULES

  Bio::EnsEMBL::ArchiveStableId
  Bio::EnsEMBL::DBSQL::ArchiveStableIdAdaptor
  Bio::EnsEMBL::StableIdHistoryTree

=cut

package Bio::EnsEMBL::StableIdEvent;

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


=head2 get_attribute

  Arg[1]      : String $type - determines whether to get attribute from 'old'
                or 'new' ArchiveStableId
  Arg[2]      : String $attr - ArchiveStableId attribute to fetch
  Example     : my $old_stable_id = $event->get_attribute('old', 'stable_id');
  Description : Accessor to attributes of the ArchiveStableIds attached to this
                event. Convenience method that does the check for undef old
                and/or new ArchiveStableId for you.
  Return type : same as respective method in Bio::EnsEMBL::ArchiveStableId, or
                undef
  Exceptions  : thrown on wrong arguments
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub get_attribute {
  my ($self, $type, $attr) = @_;

  throw("First argument passed to this function has to be 'old' or 'new'.")
    unless ($type eq 'old' or $type eq 'new');

  my %allowed_attribs = map { $_ => 1 }
    qw(stable_id version db_name release assembly);

  throw("Attribute $attr not allowed.") unless $allowed_attribs{$attr};

  my $call = $type.'_ArchiveStableId';

  if (my $id = $self->$call) {
    return $id->$attr;
  } else {
    return undef;
  }
}


=head2 ident_string

  Example     : print $event->ident_string, "\n";
  Description : Returns a string that can be used to identify your StableIdEvent.
                Useful in debug warnings.
  Return type : String
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub ident_string {
  my $self = shift;

  my $old_id = $self->old_ArchiveStableId;
  my $new_id = $self->new_ArchiveStableId;

  my $str;

  if ($old_id) {
    $str = $old_id->stable_id.'.'.$old_id->version.' ('.
      $old_id->release.')';
  } else {
    $str = 'null';
  }

  $str .= ' -> ';

  if ($new_id) {
    $str .= $new_id->stable_id.'.'.$new_id->version.' ('.
      $new_id->release.')';
  } else {
    $str .= 'null';
  }

  $str .= ' ['.$self->score.']';
  
  return $str;
}


1;

