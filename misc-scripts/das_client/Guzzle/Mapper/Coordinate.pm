
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

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Guzzle::Mapper::Coordinate;
use vars qw(@ISA);
use strict;

sub new {
  my($class, $id, $start, $end, $strand) = @_;

  return bless { 'id' => $id,
		 'start' => $start,
		 'end'   => $end,
		 'strand' => $strand }, $class;
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
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'start'} = $value;
    }
    return $obj->{'start'};

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
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'end'} = $value;
    }
    return $obj->{'end'};

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
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'strand'} = $value;
    }
    return $obj->{'strand'};

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
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'id'} = $value;
    }
    return $self->{'id'};

}

1;
