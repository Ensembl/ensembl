# EnsEMBL module for MarkerSynonym
# Copyright EMBL-EBI/Sanger center 2002
#
#
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Map::MarkerSynonym

=head1 SYNOPSIS


=head1 DESCRIPTION

Represents an alias for a marker in the EnsEMBL database. 

=cut

package Bio::EnsEMBL::Map::MarkerSynonym;

use strict;
use vars qw(@ISA);


=head2 new

  Arg [1]    : (optional) int $dbID
  Arg [2]    : (optional) string $source
  Arg [3]    : (optional) string $name
  Example    : $synonym = Bio::EnsEMBL::Map::MarkerSynonym->new(12,$src,$name);
  Description: Creates a new MarkerSynonym 
  Returntype : Bio::EnsEMBL::Map::MarkerSynonym
  Exceptions : non
  Caller     : general
  Status     : stable

=cut

sub new {
  my ($caller, $dbID, $source, $name) = @_;

  my $class = ref($caller) || $caller;

  return bless( {'dbID'   => $dbID,
		 'source' => $source,
		 'name' => $name}, $class );
}


=head2 dbID

  Arg [1]    : (optional) int $dbID
  Example    : $mid = $marker_synonym->dbID;
  Description: Getter/Setter for the internal id of this synonym
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : stable

=cut

sub dbID {
  my $self = shift;

  if(@_) {
    $self->{'dbID'} = shift;
  }

  return $self->{'dbID'};
}


=head2 source

  Arg [1]    : (optional) string $source
  Example    : $source = $marker_synonym->source;
  Description: Getter/Setter for the source of this marker synonym
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : stable

=cut

sub source {
  my $self = shift;

  if(@_) {
    $self->{'source'} = shift;
  }

  return $self->{'source'};
}


=head2 name

  Arg [1]    : (optional) string $name
  Example    : $name = $marker_synonym->name;
  Description: Getter/Setter for the name/identifier of this synonym
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : stable

=cut

sub name {
  my $self = shift;

  if(@_) {
    $self->{'name'} = shift;
  }
  
  return $self->{'name'}
}

1;

