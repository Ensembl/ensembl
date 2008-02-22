#
# Ensembl module for Bio::EnsEMBL::Mapper::Gap
#
# Written by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Mapper::Gap

=head1 SYNOPSIS

=head1 DESCRIPTION

Representation of a gap in a sequence; returned from
Mapper.pm when the target region is in a gap.

=head1 AUTHOR - Ewan Birney

This modules is part of the Ensembl project http://www.ensembl.org

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 METHODS

=cut

package Bio::EnsEMBL::Mapper::Gap;

use strict;

=head2 new

  Arg [1]    : int $start
  Arg [2]    : int $end
  Example    : $gap = Bio::EnsEMBL::Mapper::Gap($start, $end);
  Description: Creates a new Gap object.
  Returntype : Bio::EnsEMBL::Mapper::Gap
  Exceptions : none
  Caller     : Bio::EnsEMBL::Mapper
  Status     : Stable

=cut

sub new {
  my ( $proto, $start, $end ) = @_;

  my $class = ref($proto) || $proto;

  return bless( { 'start' => $start, 'end' => $end }, $class );
}

=head2 start

  Arg [1]    : (optional) int $start
               start coordinate of gap region
  Example    : $start = $gap->start();
  Description: Getter/Setter for the start attribute
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub start {
  my ( $self, $value ) = @_;

  if ( defined($value) ) {
    $self->{'start'} = $value;
  }

  return $self->{'start'};
}

=head2 end

  Arg [1]    : (optional) int $newval
               The new value to set the end coordinate to
  Example    : $end = $gap->end()
  Description: Getter/Setter for the end coordinate of the gap region
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub end {
  my ( $self, $value ) = @_;

  if ( defined($value) ) {
    $self->{'end'} = $value;
  }

  return $self->{'end'};
}

=head2 length

  Arg [1]    : none
  Example    : $len = $gap->length();
  Description: Getter for the length of this gap region
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub length {
  my ($self) = @_;

  return $self->{'end'} - $self->{'start'} + 1;
}

1;
