
#
# BioPerl module for Bio::EnsEMBL::EMBLLOAD::Obj
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::EMBLLOAD::Obj

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

package Bio::EnsEMBL::EMBLLOAD::Obj;
use vars qw(@ISA);
use strict;
use Bio::Root::Object;
@ISA = qw(Bio::Root::Object);
use Bio::AnnSeqIO;
use Bio::EnsEMBL::Gene;
use  Bio::EnsEMBL::EMBLLOAD::Clone;
use  Bio::EnsEMBL::EMBLLOAD::Contig;

sub _initialize {
    my($self,@args) = @_;
    
    my $make = $self->SUPER::_initialize;    
    my ($file,$format)=$self->_rearrange([qw(FILE FORMAT)],@args);        
# my $file='data/newentry';
# my$format='EMBL';  

   $file   || $self->throw("no file");
    $format || $self->throw("no format");

      
    my $stream= Bio::AnnSeqIO->new( -format => $format, -file => $file);
    my $annseq=$stream->next_annseq();
    $self->_get_AnnSeq($annseq);

    return $make;     
}


=head2 get_Clone

 Title   : get_Clone
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut



sub get_Clone {

    my ($self,$value)=@_;
    my $clone = Bio::EnsEMBL::EMBLLOAD::Clone->new($self->_get_AnnSeq);   
    return $clone;
}




=head2 get_Contig

 Title   : get_Contig
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Contig {

   my ($self) = @_;
   my $contig = Bio::EnsEMBL::EMBLLOAD::Contig->new($self->_get_AnnSeq);    
   return $contig;  
}





=head2 get_Gene

 Title   : get_Gene
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut


sub get_Gene {
   my ($self,$value) = @_;

   my $gene;
   my $seen;
   
   if (defined $value){
       my $contig=$self->get_Contig;
       my @genes=$contig->get_all_Genes;
       foreach my $gene(@genes){
	   if ($value eq $gene->id){$seen=1;}
	   else {$seen=0;}}}
   if ($seen==1){
       $gene = Bio::EnsEMBL::Gene->new();   
       $gene->id($value);}   
   else{$self->throw("can't get a gene without valid id");}
   
   return $gene;
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
   my $id=$self->_get_AnnSeq->seq->id;
   return $id;
   
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

=head2 donor_locator

 Title   : get_donor_locator
 Usage   : $obj->get_donor_locator; 
 Function: Reads the meta table of the database to get the donor_database_locator
 Example : get_donor_locator
 Returns : locator string
 Args    : none


=cut

sub get_donor_locator{
    my ($self) = @_;
    
    $self->throw("Not implemented in the object!");
}
=head2 get_Ghosts_by_deleted
    
 Title   : get_Ghosts_by_deleted
 Usage   : $obj->get_Ghosts_by_deleted ($recipient_last_update, $recipient_now_offset)
 Function: Gets all the ghosts for objects that have been deleted (i.e.permanently from 
	   the donor db) between the current time - offset time given by
           the recipient database and the last update time stored in its meta table 
 Example : $obj->get_Ghosts_by_deleted (973036800,973090800)
 Returns : ghost objects
 Args    : $recipient_last_update, $recipient_now_offset

=cut

sub get_updated_Ghosts{
    my ($self) = @_;
    
    $self->throw("Not implemented in the object!");
}

=head2 get_Ghost
    
 Title   : get_Ghost
 Usage   : $obj->get_Ghost ($ghost_id,$ghost_version,$ghost_obj_type)
 Function: Gets a ghost by id, version,obj_type  
 Example : $obj->get_Ghost ('test','1','transcript')
 Returns : ghost objects
 Args    : ghost id, version and object type

=cut

sub get_Ghost{
    my ($self) = @_;
    
    $self->throw("Not implemented in the object!");
}

=head2 write_Ghost
    
 Title   : write_Ghost
 Usage   : $obj->write_Ghost ($ghost)
 Function: Writes a ghost to the database  
 Example : $obj->write_Ghost ($ghost)
 Returns : 
 Args    : ghost object

=cut

sub write_Ghost{
    my ($self) = @_;
    
    $self->throw("Not implemented in the object!");
}

=head2 archive_Gene
    
 Title   : archive_Gene
 Usage   : $obj->archive_gene($gene,$clone,$arcdb)
 Function: Deletes a gene and all its transcripts and exons, 
           and archives partial info in the archive db passed on.
 Example : 
 Returns : nothing
 Args    : $gene, $arcdb (archive database object)


=cut

sub archive_Gene {
    my ($self) = @_;
 
    $self->throw("Not implemented in the object!");
}



=head2 replace_last_update
    
 Title   : replace_last_update
 Usage   : $obj->replace_last_update
 Function: Replaces the time in the last update field of the meta table with the current time
 Example : 
 Returns : nothing
 Args    : 

=cut

sub replace_last_update {
    my ($self) = @_;
 
    $self->throw("Not implemented in the object!");
}
 




sub _get_AnnSeq {

    my ($obj,$value) = @_;
    if( defined $value) {$obj->{'annseq'} = $value;}
    return $obj->{'annseq'};    
}









1;




















