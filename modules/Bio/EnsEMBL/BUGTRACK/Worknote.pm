
#
# BioPerl module for Worknote
#
# Cared for by Elia Stupka <elia@ebi.ac.uk>
#
# Copyright Elia Stupka
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Worknote - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package EnsEMBL::Worknote;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;


@ISA = qw(Bio::Root::Object);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;

# set stuff in self from @args
 return $make; # success - we hope!
}

=head2 author

 Title   : author
 Usage   : $obj->author($newval)
 Function: 
 Returns : value of author
 Args    : newvalue (optional)


=cut

sub author{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'author'} = $value;
    }
    return $obj->{'author'};
}

=head2 date

 Title   : date
 Usage   : $obj->date($newval)
 Function: 
 Returns : value of date
 Args    : newvalue (optional)


=cut

sub date{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'date'} = $value;
    }
    return $obj->{'date'};
}

=head2 note

 Title   : note
 Usage   : $obj->note($newval)
 Function: 
 Returns : value of note
 Args    : newvalue (optional)


=cut

sub note{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'note'} = $value;
    }
    return $obj->{'note'};
}
