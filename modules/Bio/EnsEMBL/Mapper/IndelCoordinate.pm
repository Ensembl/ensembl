#
# Ensembl module for Bio::EnsEMBL::Mapper::IndelCoordinate
#
# Written by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Mapper::IndelCoordinate

=head1 SYNOPSIS

=head1 DESCRIPTION

Representation of a indel in a sequence; returned from
Mapper.pm when the target region is in a deletion.

=head1 AUTHOR - Ewan Birney

This modules is part of the Ensembl project http://www.ensembl.org

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 METHODS

=cut

package Bio::EnsEMBL::Mapper::IndelCoordinate;

use Bio::EnsEMBL::Mapper::Gap;
use Bio::EnsEMBL::Mapper::Coordinate;

use vars qw(@ISA);
use strict;

@ISA = qw(Bio::EnsEMBL::Mapper::Coordinate Bio::EnsEMBL::Mapper::Gap);


=head2 new

  Arg [1]    : Bio::EnsEMBL::Mapper::Gap $gap
  Arg [2]    : Bio::EnsEMBL::Mapper::Coordinate $coordinate
  Example    : $indelCoord = Bio::EnsEMBL::Mapper::IndelCoordinate($gap, $coordinate);
  Description: Creates a new IndelCoordinate object.
  Returntype : Bio::EnsEMBL::Mapper::IndelCoordinate
  Exceptions : none
  Caller     : Bio::EnsEMBL::Mapper

=cut

sub new {
  my ( $proto, $gap, $coordinate ) = @_;

  my $class = ref($proto) || $proto;

  return
    bless( { 'start'        => $coordinate->start(),
             'end'          => $coordinate->end(),
             'strand'       => $coordinate->strand(),
             'id'           => $coordinate->id(),
             'coord_system' => $coordinate->coord_system(),
             'gap_start'    => $gap->start(),
             'gap_end'      => $gap->end()
           },
           $class );
}

=head2 gap_start

  Arg[1]      : (optional) int $gap_start
  Example     : $gap_start = $ic->gap_start()
  Description : Getter/Setter for the start of the Gap region
  ReturnType  : int
  Exceptions  : none
  Caller      : general

=cut

sub gap_start {
  my ( $self, $value ) = @_;

  if ( defined($value) ) {
    $self->{'gap_start'} = $value;
  }

  return $self->{'gap_start'};
}

=head2 gap_end

  Arg[1]      : (optional) int $gap_end
  Example     : $gap_end = $ic->gap_end()
  Description : Getter/Setter for the end of the Gap region
  ReturnType  : int
  Exceptions  : none
  Caller      : general

=cut

sub gap_end {
  my ( $self, $value ) = @_;

  if ( defined($value) ) {
    $self->{'gap_end'} = $value;
  }

  return $self->{'gap_end'};
}

=head2 gap_length

  Args        : None
  Example     : $gap_length = $ic->gap_length()
  Description : Getter for the length of the Gap region
  ReturnType  : int
  Exceptions  : none
  Caller      : general

=cut

sub gap_length {
  my ($self) = @_;

  return $self->{'gap_end'} - $self->{'gap_start'} + 1;
}

1;
