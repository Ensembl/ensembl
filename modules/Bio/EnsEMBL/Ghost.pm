
#
# BioPerl module for Ghost Object
#
# Cared for by Elia Stupka <elia@ebi.ac.uk>
#
# Copyright Elia Stupka
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Ghost - Object for Ghosts of deleted objects

=head1 SYNOPSIS

Ghost objects. 

=head1 DESCRIPTION

These are used by donor databases to pass on information about objects that have been 
permanently deleted (and archived) to recipient databases. They allow recipient databases 
to know which objects have been deleted in the donor database. These objects are 
stored in a separate table within each major database (not its archive DB), and store 
the id, type, and time of deletion for deleted objects. This is a separate concept from the 
archive db, which will hold store information about the content of the deleted objects.

=head1 CONTACT

Elia Stupka
European Bioinformatics Institute

e-mail: elia@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods 
are usually preceded with a _

=cut


# Let the code begin...
package Bio::EnsEMBL::Ghost;
use vars qw(@ISA);
use strict;

use Bio::Root::Object;

@ISA = qw(Bio::Root::Object);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize(@args);
  $self->{'_pointer_array'} = [];
# set stuff in self from @args
  return $make; # success - we hope!
}

=head2 id

 Title   : id
 Usage   : $obj->id($newval)
 Function: stores the id of the deleted object
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

=head2 version

 Title   : version
 Usage   : $obj->version($newval)
 Function: stores the version of the deleted object
 Returns : value of version
 Args    : newvalue (optional)


=cut

sub version{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'version'} = $value;
    }
    return $obj->{'version'};

}

=head2 obj_type

 Title   : obj_type
 Usage   : $obj->obj_type($newval)
 Function: stores the object type of the deleted object
 Returns : value of obj_type
 Args    : newvalue (optional)


=cut

sub obj_type{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'obj_type'} = $value;
    }
    return $obj->{'obj_type'};

}

=head2 deleted

 Title   : deleted
 Usage   : $obj->deleted($newval)
 Function: stores the time of deletion of the deleted object
 Returns : value of deleted
 Args    : newvalue (optional)


=cut

sub deleted{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'deleted_time'} = $value;
    }
    return $obj->{'deleted_time'};

}

=head2 add_pointer_id

 Title   : add_pointer_id
 Usage   : $ghost->add_pointer_id('ENSG0000000012');
 Function: Adds the id of a forward reference, for the new
           object (if present) of this deleted object. There
           is no ability to put in new classes of objects
 Returns : nothing
 Args    : an id


=cut

sub add_pointer_id{
   my ($self,$value) = @_;

   push(@{$self->{'_pointer_array'}},$value);
}

=head2 each_pointer_id

 Title   : each_pointer_id
 Usage   : foreach $id ( $ghost->each_pointer_id() )
 Function: returns all the forward references for the new
           object of this deleted object.
 Example :
 Returns : an array of strings
 Args    :


=cut

sub each_pointer_id{
   my ($self) = @_;

   return @{$self->{'_pointer_array'}};
}


=head2 _stored

 Title   : _stored
 Usage   : $obj->_stored($newval)
 Function: Internal method, should not really be needed
           stores the time of storage of the deleted object
 Returns : value of _stored
 Args    : newvalue (optional)


=cut

sub _stored{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_stored'} = $value;
    }
    return $obj->{'_stored'};
}
