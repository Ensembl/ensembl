
#
# BioPerl module for DB/ContigI.pm
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DB::ContigI.pm - Abstract Interface for Contig

=head1 SYNOPSIS

This is the abstract definition of a Contig, along with 'decorator'
functions

=head1 DESCRIPTION

Describe the object here

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::EMBLLOAD::Contig;
use vars qw(@ISA);
use strict;
use Bio::Root::Object;
use Bio::AnnSeq;
@ISA = qw(Bio::Root::Object Bio::EnsEMBL::DB::ContigI);
use Bio::EnsEMBL::EMBLLOAD::Obj;
use Bio::EnsEMBL::Analysis::Analysis;
use Bio::EnsEMBL::SeqFeature;

sub _initialize {
    my($self,@args) = @_;
    
    my ($annseq)=$self->_rearrange([qw(ANNSEQ)],@args);
    my $make = $self->SUPER::_initialize; 
    $self->_get_AnnSeq($annseq); 
    $self->id;
    return $make; 
    
}



=head2 id

 Title   : id
 Usage   : $obj->id($newval)
 Function: 
 Example : 
 Returns : value of id
 Args    : newvalue (optional)


=cut

sub id {

   my ($self) = @_;
   my  $id=$self->_get_AnnSeq->seq->id . "00001"; 
   return $id;

}



=head2 seq

 Title   : seq
 Usage   : $seq = $contig->seq();
 Function: Gets a Bio::Seq object out from the contig
 Example :
 Returns : Bio::Seq object
 Args    :


=cut

sub seq {
   my ($self) = @_;
   my $seq=$self->_get_AnnSeq->seq;
   return $seq;

}




=head2 get_all_SeqFeatures

 Title   : get_all_SeqFeatures
 Usage   : foreach my $sf ( $contig->get_all_SeqFeatures ) 
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_SeqFeatures {
    my ($self) = @_;
    my  @features=$self->_get_AnnSeq->all_SeqFeatures;	
    my @ensembl_features;
    
    foreach my $feature(@features){
	
	# Ewan to explain why do I have to copy one object to another

	my $analysis = new Bio::EnsEMBL::Analysis::Analysis(-db              => 'EMBL',
							    -db_version      => 'NULL',
							    -program         => 'NULL',
							    -program_version => 'NULL',
							    -gff_source      => 'NULL',
							    -gff_feature     => 'EMBL ann',
							    );
	
	my $ensembl_feature=new Bio::EnsEMBL::SeqFeature(-seqname => $self->id,
							 -start   => $feature->start,
							 -end     => $feature->end,
							 -strand  => $feature->strand,
							 -frame   => $feature->frame,
							 -source_tag  => $feature->source_tag,
							 -primary_tag => $feature->primary_tag,
							 -analysis => $analysis,
							 -score => $feature->score
							 );	
	push @ensembl_features,$ensembl_feature;
    }
    
    return @ensembl_features;  
}



=head2 get_all_SimilarityFeatures

 Title   : get_all_SimilarityFeatures
 Usage   : foreach my $sf ( $contig->get_all_SimilarityFeatures ) 
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_SimilarityFeatures{
   my ($self) = @_;
   my @sim_features;
   my @features=$self->get_all_SeqFeatures;
   foreach my $feature (@features){
       if ($feature->analysis->gff_feature eq 'similarity'){
	   push @sim_features,$feature;}}

       return @sim_features;
}




=head2 get_all_Genes

 Title   : get_all_Genes
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut



sub get_all_Genes {
	
    my ($self)=@_;
    my @genes;
    my @exons;
    my $exon_counter;    
    foreach my $ft($self->_get_AnnSeq->all_SeqFeatures){
	if($ft->primary_tag eq 'CDS'){	    
	    my $exon = Bio::EnsEMBL::Exon->new($ft->start,$ft->end,$ft->strand);
	    $exon->phase("1");
	    $exon->end_phase("1");
	    $exon_counter++;
	    $exon->id($exon_counter);	
	    $exon->contig_id($self->id);
	    $exon->version("1");
	    $exon->created("2000");
	    $exon->modified("2000");
	    $exon->attach_seq($self->_get_AnnSeq->seq);
   
	    push @exons,$exon;
	}	
    }
    
    unless (scalar @exons ==0){
	my $transcript = Bio::EnsEMBL::Transcript->new(@exons);
	$transcript->id("transcript_id");	
	my $translation=Bio::EnsEMBL::Translation->new();
	$translation->id ("some_id");
	$translation->version (2);
	$translation->start (55);
	$translation->start_exon_id (1);
	$translation->end (55);
	$translation->end_exon_id (3);

	$transcript->version(1);	
	#$transcript->gene("new_gene_id");	
	$transcript->translation($translation);	
	my $gene = Bio::EnsEMBL::Gene->new();   
	my $gene_id=$self->_get_AnnSeq->seq->id;
	$gene_id="EMBLG" . "0000" . $gene_id;
	$gene->id($gene_id);
	$gene->add_Transcript($transcript);
	
	push @genes,$gene;
    }
    return @genes;
}






=head2 length

 Title   : length
 Usage   : 
 Function: Provides the length of the contig
 Example :
 Returns : 
 Args    :


=cut

sub length {
   my ($self,@args) = @_;
   my $length=$self->_get_AnnSeq->seq->seq_len;
   return $length;

}






sub _get_AnnSeq {
    my ($self,$value) = @_;
    if (defined $value){$self->{'annseq'}=$value;}
    return $self->{'annseq'};
}







=head2 offset

 Title   : offset
 Usage   : $offset = $contig->offset()
 Function: Provides the offset of the contig in the clone
         : somehow. 1 means it is the first contig
 Example :
 Returns : 
 Args    :


=cut



sub offset{
   my ($self,@args) = @_;

   my $offset=2;
   #$self->throw("Object did not provide the offset method on Contig interface!");
   return $offset;
}



=head2 created

 Title   : created
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub created{
   my ($self) = @_;
   #$self->throw("Class [$self] has not implemented the created method");
   my $created=4;
   return $created;

}




=head2 seq_date

 Title   : seq_date
 Usage   : $contig->seq_date()
 Function: Gives the unix time value of the dna table created datetime field, which indicates
           the original time of the dna sequence data
 Example : $contig->seq_date()
 Returns : unix time
 Args    : none


=cut

sub seq_date{
    my ($self) = @_;

    my $seq_date=22;
#    $self->throw("Object did not provide the seq_date method on Contig interface!");
    return $seq_date;
}



=head2 orientation

 Title   : orientation
 Usage   : 
 Function: Provides the orientation of the contig in the clone.
 Example :
 Returns : 
 Args    :


=cut

sub orientation{
   my ($self,@args) = @_;

   #$self->throw("Object did not provide the orientation method on Contig interface!");

   my $orientation=1;
   return $orientation;

}

=head2 order

 Title   : order
 Usage   : $obj->order($newval)
 Function: 
 Returns : value of order
 Args    : newvalue (optional)


=cut

sub order{
    my ($self,@args) = @_;

   # $self->throw("Object did not provide the order method on Contig interface!");
    my $order=1;
    return $order;

}





1;
