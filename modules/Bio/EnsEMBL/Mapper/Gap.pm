

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

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Mapper::Gap;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::EnsEMBL::Root

use Bio::EnsEMBL::Root;

@ISA = qw(Bio::EnsEMBL::Root);

# new() is written here 

sub new {
  my($class,$start, $end) = @_;

  return bless { 'start' => $start,
		 'end'   => $end }, $class;

}

=head2 start

  Arg  1      int $start
	      start coordinate of gap region
  Function    accessor method
  Returntype  int
  Exceptions  none
  Caller      Bio::EnsEMBL::Mapper::Gap

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
	      end coordinate of gap region
  Function    accessor method
  Returntype  int
  Exceptions  none
  Caller      Bio::EnsEMBL::Mapper::Gap

=cut

sub end{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'end'} = $value;
    }
    return $self->{'end'};

}

sub length {
  my $self = shift;
  return $self->{'end'} - $self->{'start'} + 1;
}



1;
