
#
# Ensembl module for Bio::EnsEMBL::MapperCoordinate
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::MapperCoordinate - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=head1 CONTACT

Ensembl - ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::MapperCoordinate;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::RootI

use Bio::Root::RootI;


@ISA = qw(Bio::Root::RootI);

sub new {
  my($class,@args) = @_;

  my $self = {};
  bless $self,$class;
  
  return $self;
}

=head2 start

 Title   : start
 Usage   : $obj->start($newval)
 Function: 
 Returns : value of start
 Args    : newvalue (optional)


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

 Title   : end
 Usage   : $obj->end($newval)
 Function: 
 Returns : value of end
 Args    : newvalue (optional)


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

 Title   : strand
 Usage   : $obj->strand($newval)
 Function: 
 Returns : value of strand
 Args    : newvalue (optional)


=cut

sub strand{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'strand'} = $value;
    }
    return $obj->{'strand'};

}

=head2 abutts

 Title   : abutts
 Usage   : $obj->abutts($newval)
 Function: 
 Returns : value of abutts
 Args    : newvalue (optional)


=cut

sub abutts{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'abutts'} = $value;
    }
    return $obj->{'abutts'};

}

=head2 has_orig_start

 Title   : has_orig_start
 Usage   : $obj->has_orig_start($newval)
 Function: 
 Returns : value of has_orig_start
 Args    : newvalue (optional)


=cut

sub has_orig_start{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'has_orig_start'} = $value;
    }
    return $obj->{'has_orig_start'};

}

=head2 has_orig_end

 Title   : has_orig_end
 Usage   : $obj->has_orig_end($newval)
 Function: 
 Returns : value of has_orig_end
 Args    : newvalue (optional)


=cut

sub has_orig_end{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'has_orig_end'} = $value;
    }
    return $obj->{'has_orig_end'};

}

1;
