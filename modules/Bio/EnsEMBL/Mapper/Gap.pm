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

Bio::EnsEMBL::Mapper::Gap

=head1 SYNOPSIS

=head1 DESCRIPTION

Representation of a gap in a sequence; returned from Mapper.pm when the
target region is in a gap.

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
  my ( $proto, $start, $end, $rank ) = @_;

  my $class = ref($proto) || $proto;

  return bless( { 'start' => $start, 'end' => $end, 'rank' => $rank  || 0 }, $class );
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

sub rank {
  my ( $self, $value ) = @_;

  if ( defined($value) ) {
    $self->{'rank'} = $value;
  }

  return $self->{'rank'};
}

1;
