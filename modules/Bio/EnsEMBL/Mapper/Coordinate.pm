
#
# Ensembl module for Bio::EnsEMBL::MapperCoordinate
#
# Written by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Mapper::Coordinate

=head1 SYNOPSIS

=head1 DESCRIPTION

Representation of a mapped region in a sequence; returned from
Mapper.pm when the target region maps on to valid sequence.

=head1 CONTACT

This module is part of the Ensembl project http://www.ensembl.org

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 METHODS

=cut

package Bio::EnsEMBL::Mapper::Coordinate;
use vars qw(@ISA);
use strict;

sub new {
  my($class, $id, $start, $end, $strand, $coord_system) = @_;

  return bless { 'id' => $id,
                 'start' => $start,
                 'end'   => $end,
                 'strand' => $strand,
                 'coord_system' => $coord_system}, $class;
}


=head2 start

  Arg  1      int $start
              start coordinate of object in mapped region
  Function    accessor method
  Returntype  int
  Exceptions  none
  Caller      Bio::EnsEMBL::Mapper::Coordinate

=cut

sub start{
  my $self = shift;
  $self->{'start'} = shift if(@_);
  return $self->{'start'};
}


=head2 end

  Arg  1      int $end
              end coordinate of object in mapped region
  Function    accessor method
  Returntype  int
  Exceptions  none
  Caller      Bio::EnsEMBL::Mapper::Coordinate

=cut

sub end{
  my $self = shift;
  $self->{'end'} = shift if(@_);
  return $self->{'end'};
}


=head2 strand

  Arg  1      int $strand
              strand of object in mapped region
  Function    accessor method
  Returntype  int
  Exceptions  none
  Caller      Bio::EnsEMBL::Mapper::Coordinate

=cut

sub strand{
  my $self = shift;
  $self->{'strand'} = shift if(@_);
  return $self->{'strand'};
}


=head2 id

  Arg  1      char|int $id
              id of object in mapped region
	      e.g. RawContig ID/chromosome name
  Function    accessor method
  Returntype  char|int
  Exceptions  none
  Caller      Bio::EnsEMBL::Mapper::Coordinate

=cut

sub id{
  my $self = shift;
  $self->{'id'} = shift if(@_);
  return $self->{'id'};
}


sub coord_system {
  my $self = shift;
  $self->{'coord_system'} = shift if(@_);
  return $self->{'coord_system'};
}

sub length {
  my $self = shift;
  return $self->{'end'} - $self->{'start'} + 1;
}

1;
