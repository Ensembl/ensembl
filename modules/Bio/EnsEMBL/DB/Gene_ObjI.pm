#
# EnsEMBL module for Bio::EnsEMBL::DBSQL::Gene_Obj
#
# Cared for by Elia Stupka <elia@ebi.ac.uk>
#
# Copyright Elia Stupka
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DB::Gene_ObjI - MySQL database adapter interface for EnsEMBL genes, transcripts
exons, etc.

=head1 SYNOPSIS


=head1 DESCRIPTION

This is one of the interfaces contained in Bio:EnsEMBL::DB::ObjI, dealing with
Gene methods, such as writing and getting genes, transcripts, translations, and exons.

=head1 CONTACT

Elia Stupka: elia@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are 
usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::DBSQL::Gene_ObjI;

use vars qw(@ISA);
use strict;


=head2 get_all_my_cloneid

 Title   : get_all_my_cloneid
 Usage   : @each_cloneid = $gene_obj->get_all_my_cloneid($geneid);
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_my_cloneid{
    my ($self) = @_;

   $self->throw("Not implemented in the object!");

}

=head2 delete

 Title   : delete
 Usage   : $Gene_Obj->delete_Gene($gene_id)
 Function: deletes a gene from the database, i.e. exons, transcripts, translations
 Example : $geneobj->delete_Gene('ENSG00000019482')
 Returns : nothing
 Args    : gene id

=cut

sub delete{
    my ($self) = @_;

   $self->throw("Not implemented in the object!");

}   

=head2 delete_Exon

 Title   : delete_Exon
 Usage   : $obj->delete_Exon($exon_id)
 Function: Deletes exon, including exon_transcript rows
 Example : $obj->delete_Exon(ENSE000034)
 Returns : nothing
 Args    : $exon_id

=cut

sub delete_Exon{
    my ($self) = @_;

   $self->throw("Not implemented in the object!");
}

=head2 delete_Supporting_Evidence

 Title   : delete_Supporting_Evidence
 Usage   : $obj->delete_Supporting_Evidence($exon_id)
 Function: Deletes exon\'s supporting evidence entries
 Example : $obj->delete_Supporting_Evidence(ENSE000034)
 Returns : nothing
 Args    : $exon_id


=cut

sub delete_Supporting_Evidence {
    my ($self) = @_;

   $self->throw("Not implemented in the object!");

}

=head2 get

 Title   : get
 Usage   : $geneobj->get($geneid, $supporting)
 Function: gets one gene out of the db with or without supporting evidence
 Example : $obj->get('ENSG00000009151','evidence')
 Returns : gene object (with transcripts, exons and supp.evidence if wanted)
 Args    : gene id and supporting tag (if latter not specified, assumes without
	   Note that it is much faster to get genes without supp.evidence!

=cut

sub get {
    my ($self) = @_;

   $self->throw("Not implemented in the object!");

}

=head2 get_all_Gene_id

 Title   : get_all_Gene_id
 Usage   : $geneobj->get_all_Gene_id
 Function: Gets an array of ids for all genes in the current db
 Example : $geneobj->get_all_Gene_id
 Returns : array of ids
 Args    : none

=cut

sub get_all_Gene_id{
    my ($self) = @_;

   $self->throw("Not implemented in the object!");

}

=head2 get_array_supporting

    Title   : get_Gene_array_supporting
    Usage   : $obj->get_Gene_array_supporting($supporting,@geneid)
    Function: Gets an array of genes, with transcripts and exons. If $supporting
           equal to 'evidence' the supporting evidence for each exon is also read
    from the supporting evidence table
    Example : $obj->get_Gene_array_supporting ('evidence',@geneid)
    Returns : an array of gene objects
    Args    : 'evidence' and gene id array

    
=cut
    
sub get_array_supporting {
    my ($self) = @_;

   $self->throw("Not implemented in the object!");

}

=head2 get_Exon

 Title   : get_Exon
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Exon{
    my ($self) = @_;

   $self->throw("Not implemented in the object!");
}

=head2 get_supporting_evidence

 Title   : get_supporting_evidence
 Usage   : $obj->get_supporting_evidence
 Function: Writes supporting evidence features to the database
 Example :
 Returns : nothing
 Args    : array of exon objects, needed to know which exon to attach the evidence to


=cut

sub get_supporting_evidence {
    my ($self) = @_;

   $self->throw("Not implemented in the object!");
}

=head2 get_Transcript
    
 Title   : get_Transcript
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut
    
sub get_Transcript{
    my ($self) = @_;

   $self->throw("Not implemented in the object!");

}

=head2 get_Translation

 Title   : get_Translation
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Translation{
    my ($self) = @_;

   $self->throw("Not implemented in the object!");

}

=head2 write

 Title   : write
 Usage   : $Gene_obj->write_Gene($gene)
 Function: writes a particular gene into the database
 Example :
 Returns : nothing
 Args    : $gene object


=cut

sub write{
    my ($self) = @_;

   $self->throw("Not implemented in the object!");

}

=head2 write_Exon

 Title   : write_Exon
 Usage   : $obj->write_Exon($exon)
 Function: writes a particular exon into the database
 Example :
 Returns : 
 Args    :


=cut

sub write_Exon {
    my ($self) = @_;

   $self->throw("Not implemented in the object!");
}

=head2 write_supporting_evidence

 Title   : write_supporting_evidence
 Usage   : $obj->write_supporting_evidence
 Function: Writes supporting evidence features to the database
 Example :
 Returns : nothing
 Args    : None


=cut

sub write_supporting_evidence {
    my ($self) = @_;

   $self->throw("Not implemented in the object!");
}

=head2 write_Transcript

 Title   : write_Transcript
 Usage   : $obj->write_Transcript($trans,$gene)
 Function: writes a particular transcript *but not the exons* into
           the database
 Example :
 Returns : 
 Args    :


=cut

sub write_Transcript{
    my ($self) = @_;

   $self->throw("Not implemented in the object!");

}

=head2 write_Translation

 Title   : write_Translation
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub write_Translation { 
    my ($self) = @_;

   $self->throw("Not implemented in the object!");


}

=head2 get_new_GeneID

 Title   : get_new_GeneID
 Usage   : my $id = $geneobj->get_new_GeneID
 Function: 
 Example : 
 Returns : Gets the next unused gene id from the database
 Args    : none


=cut

sub get_new_GeneID {
    my ($self) = @_;

   $self->throw("Not implemented in the object!");

}

=head2 get_new_TranscriptID

 Title   : get_new_TranscriptID
 Usage   : my $id = $geneobj->get_new_TranscriptID
 Function: 
 Example : 
 Returns : Gets the next unused transcript id from the database
 Args    : none


=cut

sub get_new_TranscriptID {

    my ($self) = @_;

   $self->throw("Not implemented in the object!");

}

=head2 get_new_ExonID

 Title   : get_new_ExonID
 Usage   : my $id = $geneobj->get_new_ExonID
 Function: 
 Example : 
 Returns : Gets the next unused exon id from the database
 Args    : none


=cut

sub get_new_ExonID {
    my ($self) = @_;

   $self->throw("Not implemented in the object!");

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


1;
