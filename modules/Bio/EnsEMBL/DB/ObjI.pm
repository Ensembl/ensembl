
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

=head2 write_Clone

 Title   : write_Clone
 Usage   : $obj->write_Clone($cloneid,$dna)
 Function: writes a Clone and its dna into the database
 Example :
 Returns : 
 Args    :


=cut

sub write_Clone {
   my ($self) = @_;

   $self->throw("Not implemented in the object!");
}

=head2 add_ExternalFeatureFactory

 Title   : add_ExternalFeatureFactory (Abstract)
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_ExternalFeatureFactory{
   my ($self) = @_;
   $self->throw("Abstract method add_ExternalFeatureFactory encountered in base class. Implementation failed to complete it")

}


1;
