
#
# BioPerl module for Bio::EnsEMBL::DB::ObjI
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DB::ObjI - Abstract Interface of Database objects generically

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

package Bio::EnsEMBL::DB::ObjI;
use vars qw($AUTOLOAD @ISA);
use strict;

=head2 get_Gene

 Title   : get_Gene
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Gene{
   my ($self) = @_;

   $self->throw("Not implemented in the object!");
}

=head2 get_Clone

 Title   : get_Clone
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Clone{
   my ($self) = @_;

   $self->throw("Not implemented in the object!");
}


=head2 get_Contig

 Title   : get_Contig
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Contig{
   my ($self) = @_;

   $self->throw("Not implemented in the object!");
}

=head2 get_all_Clone_id

 Title   : get_all_Clone_id
 Usage   : @cloneid = $obj->get_all_Clone_id
 Function: returns all the valid (live) Clone ids in the database
 Example :
 Returns : 
 Args    :


=cut

sub get_all_Clone_id{
   my ($self) = @_;
   
   $self->throw("Not implemented in the object!");
   
}

=head2 write_Gene

 Title   : write_Gene
 Usage   : $obj->write_Gene($gene)
 Function: writes a particular gene into the database
           
 Example :
 Returns : 
 Args    :


=cut

sub write_Gene{
   my ($self) = @_;

   $self->throw("Not implemented in the object!");
}

=head2 write_Contig

 Title   : write_Contig
 Usage   : $obj->write_Contig($contigid,$dna)
 Function: writes a contig and its dna into the database
 Example :
 Returns : 
 Args    :


=cut

sub write_Contig {
   my ($self) = @_;

   $self->throw("Not implemented in the object!");
}

=head2 get_updated_objects
    
 Title   : get_updated_objects
 Usage   : $obj->get_updated_objects ($recipient_last_update, $recipient_now, $recipient_offset)
 Function: Gets all the objects that have been updated (i.e.change in 
	   version number) between the current time - offset time given by
           the recipient database and the last update time stored in its meta table 
 Example : $obj->get_updated_objects (973036800,973090800)
 Returns : all the objects updated within that timespan
 Args    : $recipient_last_update, $recipient_now

=cut

sub get_updated_Objects{
    my ($self) = @_;
    
   $self->throw("Not implemented in the object!");
}

=head2 get_last_update

 Title   : get_last_update
 Usage   : $obj->get_last_update; 
 Function: Reads the meta table of the database to get the last_update time
 Example : get_last_update
 Returns : UNIX TIME of last update
 Args    : none


=cut

sub get_last_update{
    my ($self) = @_;
    
    $self->throw("Not implemented in the object!");
}

=head2 get_now_offset

 Title   : get_now_offset
 Usage   : $obj->get_now_minus_offset; 
 Function: Gets the current time from the point of view of the database, substracts the
           offset time found in the meta table and gives back unix time of now-offset
 Example : get_now_offset
 Returns : UNIX TIME of now - offset_time
 Args    : none


=cut

sub get_now_offset{
    my ($self) = @_;
    
    $self->throw("Not implemented in the object!");
} 

1;
