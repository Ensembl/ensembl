
#
# BioPerl module for Bug
#
# Cared for by Elia Stupka<elia@ebi.ac.uk>
#
# Copyright Elia Stupka
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bug - Object for the representation of BUGTRACK bugs

=head1 SYNOPSIS

Bugs relating to the ENSEMBL modules or web pages, to be stored
in the bugtrack database, and to be created through the web pages

=head1 DESCRIPTION

Bug object

=head1 CONTACT

Elia Stupka
e-mail: elia@ebi.ac.uk
phone: +44 1223 494431

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package EnsEMBL::Bug;
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
  $self->{'_worknote_array'} = [];
# set stuff in self from @args
 return $make; # success - we hope!
}

=head2 id

 Title   : id
 Usage   : $obj->id($newval)
 Function: Get/set method for the bug id
 Returns : value of id
 Args    : newvalue (optional)


=cut

sub id{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'id'} = $value;
    }
    return $obj->{'id'};

}

=head2 title

 Title   : title
 Usage   : $obj->title($newval)
 Function: 
 Returns : value of title
 Args    : newvalue (optional)


=cut

sub title{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'title'} = $value;
    }
    return $obj->{'title'};

}

=head2 add_Worknote

 Title   : add_Worknote
 Usage   : $bug->add_Worknote($note)
 Function: Adds a worknote object to the Bug
 Example : $bug->add_Worknote($note)
 Returns : nothing
 Args    : Worknote object


=cut

sub add_Worknote{
   my ($self,$note) = @_;

   if( ! $note->isa("Bio::EnsEMBL::BUGTRACK::Worknote") ) {
       $self->throw("$note is not a Bio::EnsEMBL::BUGTRACK::Worknote object!");
   }

   push(@{$self->{'_worknote_array'}},$note);
}

=head2 each_Worknote

 Title   : each_Worknote
 Usage   : foreach $trans ( $gene->each_Worknote)
 Function:
 Example :
 Returns : An array of Worknote objects
 Args    :


=cut

sub each_Worknote {
   my ($self) = @_;
   return @{$self->{'_worknote_array'}};   

}

=head2 first_Worknote

 Title   : first_Worknote
 Usage   : $bug->first_Worknote
 Function: Returns the last Worknote in the Worknote array of this bug
 Example : $bug->first_Worknote
 Returns : Worknote object
 Args    : 


=cut

sub first_Worknote{
   my ($self) = @_;
   my @temp = @{$self->{'_worknote_array'}};

   return shift @temp;
}

=head2 last_Worknote

 Title   : last_Worknote
 Usage   : $bug->last_Worknote
 Function: Returns the last Worknote in the Worknote array of this bug
 Example : $bug->last_Worknote
 Returns : Worknote object
 Args    : 


=cut

sub last_Worknote{
   my ($self) = @_;
   my @temp = @{$self->{'_worknote_array'}};

   return pop @temp;
}

=head2 type

 Title   : type
 Usage   : $obj->type($newval)
 Function: 
 Returns : value of type
 Args    : newvalue (optional)


=cut

sub type{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'type'} = $value;
    }
    return $obj->{'type'};
}
