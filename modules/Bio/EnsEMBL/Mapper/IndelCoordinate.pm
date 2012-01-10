=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::Mapper::IndelCoordinate

=head1 SYNOPSIS

=head1 DESCRIPTION

Representation of a indel in a sequence; returned from Mapper.pm when
the target region is in a deletion.

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
