
#
# Ensembl module for Bio::EnsEMBL::Mapper::Unit
#
# Written by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Mapper::Unit - One side of a map pair

=head1 SYNOPSIS

=head1 DESCRIPTION
 
Two regions mapped between different coordinate systems are
each represented by a Bio::EnsEMBL::Mapper::Unit and joined
together as a Bio::EnsEMBL::Mapper::Pair.

=head1 AUTHOR - Ewan Birney

This module is part of the Ensembl project http://www.ensembl.org

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 METHODS

=cut


package Bio::EnsEMBL::Mapper::Unit;

use strict;

sub new {
  my ( $proto, $id, $start, $end ) = @_;

  my $class = ref($proto) || $proto;

  return
    bless( { 'id' => $id, 'start' => $start, 'end' => $end }, $class );
}

=head2 id

  Arg  1      int|char $id
	      the id of the object (e.g. seq_region_name) which is mapped
  Function    accessor method
  Returntype  int|char
  Exceptions  none
  Caller      Bio::EnsEMBL::Mapper::Unit
  Status      Stable

=cut

sub id {
  my ( $self, $value ) = @_;

  if ( defined($value) ) {
    $self->{'id'} = $value;
  }

  return $self->{'id'};
}

=head2 start

  Arg  1      int $start
	      the start coordinate of the mapped
	      region which this object represents
  Function    accessor method
  Returntype  int
  Exceptions  none
  Caller      Bio::EnsEMBL::Mapper::Unit
  Status      Stable

=cut

sub start {
  my ( $self, $value ) = @_;

  if ( defined($value) ) {
    $self->{'start'} = $value;
  }

  return $self->{'start'};
}

=head2 end

  Arg  1      int $end
	      the end coordinate of the mapped
	      region which this object represents
  Function    accessor method
  Returntype  int
  Exceptions  none
  Caller      Bio::EnsEMBL::Mapper::Unit
  Status      Stable

=cut

sub end {
  my ( $self, $value ) = @_;

  if ( defined($value) ) {
    $self->{'end'} = $value;
  }

  return $self->{'end'};
}

1;
