#
# EnsEMBL module for Bio::EnsEMBL::TimDB::Update_Obj
#
# Cared for by Elia Stupka <elia@ebi.ac.uk>
#
# Copyright Elia Stupka
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::TimDB::Update_Obj - MySQL database adapter class for EnsEMBL update system

=head1 SYNOPSIS

  use Bio::EnsEMBL::TimDB::Obj;
  use Bio::EnsEMBL::TimDB::Update_Obj;

  $db = new Bio::EnsEMBL::TimDB::Obj( -user => 'root', -db => 'pog' , -host => 'caldy' , -driver => 'mysql' );
  my $update_obj=Bio::EnsEMBL::Update_Obj->new($obj);

  # Get the last update time - offset
  $update_obj->get_last_update_offset();

=head1 DESCRIPTION

This is one of the objects contained in Bio:EnsEMBL::TimDB::Obj, dealing with
the update system, such identifying last update, getting updated objects, ghosts, etc.

The Obj object represents a database that is implemented somehow (you shouldn\'t
care much as long as you can get the object). 

=head1 CONTACT

Elia Stupka: elia@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are 
usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::TimDB::Update_Obj;

use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;
use Bio::EnsEMBL::TimDB::Obj;

@ISA = qw(Bio::Root::Object);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,$db_obj) = @_;

  my $make = $self->SUPER::_initialize;
  
  $db_obj || $self->throw("Database Gene object must be passed a db obj!");
  $self->_db_obj($db_obj);

  return $make; # success - we hope!
}

=head2 get_updated_Clone_id

 Title   : get_updated_Clone_id
 Usage   : @cloneid = $obj->get_updated_Clone_id($last,$now_offset,flag)
 Function: returns all the valid (live) Clone ids in the database
 Example :
 Returns : 
 Args    : if $flag set, returns all clones regardless of invalid SV

=cut

sub get_updated_Clone_id{
    my ($self,$last,$now_offset,$fall) = @_;
    my @objs;

    # FIXME - time offset of 30mins = 30x60 seconds
    my $offset_time = 1800;
    
    $last = $last - $offset_time;
    
    # get list of updated clones
    my @clones;

    my($val,$key);

    while(($key,$val) = each %{$self->_db_obj->{'_clone_update_dbm'}}){

	my($date2)=split(',',$val);
	# make list of updatable clones
	
	if($date2 > $last && $date2 < $now_offset){
	    push(@clones,$key);
	}
    }
    
    # check validity of clones selected, then fetch
    return $self->_db_obj->_get_Clone_id($self,$fall,\@clones);
}

=head2 _db_obj

 Title   : _db_obj
 Usage   : $obj->_db_obj($newval)
 Function: 
 Example : 
 Returns : value of _db_obj
 Args    : newvalue (optional)


=cut

sub _db_obj{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_db_obj'} = $value;
    }
    return $self->{'_db_obj'};

}
