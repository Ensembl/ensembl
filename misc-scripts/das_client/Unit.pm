
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

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Unit;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::EnsEMBL::Root

# new() is written here 

sub new {
  my($class,@args) = @_;

    my $self = {};
    bless $self,$class;

# set stuff in self from @args
    return $self;
}


=head2 id

  Arg  1      int|char $id
	      the id of the object (e.g. chromosome
	      or RawContig) which is mapped
  Function    accessor method
  Returntype  int|char
  Exceptions  none
  Caller      Bio::EnsEMBL::Mapper::Unit

=cut

sub id{
   my ($self,$value) = @_;
   if( defined $value) {
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

=cut

sub start{
   my ($self,$value) = @_;
   if( defined $value) {
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

=cut

sub end{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'end'} = $value;
    }
    return $self->{'end'};

}


1;
