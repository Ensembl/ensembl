
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
use Bio::Root::RootI;

@ISA = qw(Bio::Root::RootI Bio::EnsEMBL::DB::ContigI);
use Bio::EnsEMBL::EMBLLOAD::Obj;
use Bio::EnsEMBL::Gene;


sub new {
    my($class,@args) = @_;
    
    my $self = {};
    bless $self,$class;

    my ($annseq)=$self->_rearrange([qw(ANNSEQ)],@args);

    $self->_get_Seq($annseq); 
    $self->id;
    return $self; 
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
   my  $id=$self->_get_Seq->id . ".00001"; 
   return $id;

}


=head2 internal_id

 Title   : internal_id
 Usage   : $obj->internal_id($newval)
 Function: 
 Example : 
 Returns : value of internal_id
 Args    : newvalue (optional)


=cut

sub internal_id{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'internal_id'} = $value;
    }
    return $obj->{'internal_id'};

}



=head2 seq

 Title   : seq
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub seq{
   my ($self,@args) = @_;

   return $self->primary_seq->seq;
}

=head2 primary_seq

 Title   : primary_seq
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub primary_seq{
   my ($self,@args) = @_;

   my $seq=$self->_get_Seq;
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
    return ();
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

   return ();
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

    my $exoncounter = 1;    
    my $genecounter = 1;
    
    my $id = $self->_get_Seq->id;

    my $time = time();

    foreach my $ft ( $self->_get_Seq->top_SeqFeatures ) {
	if( $ft->primary_tag eq 'CDS_span' ) {
	    #print STDERR "multi exon gene!\n";

	    my $gene       = Bio::EnsEMBL::Gene->new();
	    my $trans      = Bio::EnsEMBL::Transcript->new();
	    $gene->add_Transcript($trans);
	    $gene->id($id.".gene.".$genecounter);
	    $gene->version(1);
	    $trans->id($id.".trans.".$genecounter);
	    $trans->version(1);

	    if( $ft->has_tag('pseudo') ) {
		$gene->type('pseudo');
	    } else {
		$gene->type('standard');
	    }

	    # split seqfeature
	    my $phase = 0;
	    foreach my $sub ( $ft->sub_SeqFeature ) {
		my $exon = Bio::EnsEMBL::Exon->new();
		$exon->phase($phase);
		$exon->start($sub->start);
		$exon->end($sub->end);
		$exon->strand($sub->strand);
		$exon->contig_id($self->id);
		$exon->seqname($self->id);
		$exon->version(1);
		$exon->created($time);
		$exon->modified($time);
		$exon->id($id.".exon.".$exoncounter++);
		$trans->add_Exon($exon);
		$phase = $exon->end_phase();

	    }
	    my @exons = $trans->each_Exon;
	    my $first = shift @exons;
	    my $last;
	    if( $#exons == -1 ) {
		$last = $first;
	    } else {
		$last = pop @exons;
	    }

	    my $tranl = Bio::EnsEMBL::Translation->new();
	    $tranl->id($id.".transl.".$genecounter);
	    $tranl->start_exon_id($first->id);
	    $tranl->end_exon_id($last->id);
	    $tranl->start(1);
	    $tranl->end($last->length);
	    $tranl->version(1);
	    $trans->translation($tranl);
	    
	    $genecounter++;
	    push(@genes,$gene);
	} elsif ( $ft->primary_tag eq 'CDS' ) {
	    #print STDERR "Single exon gene!\n";

	    
	    my $gene       = Bio::EnsEMBL::Gene->new();
	    my $trans      = Bio::EnsEMBL::Transcript->new();
	    $gene->add_Transcript($trans);
	    $gene->version(1);
	    $gene->id($id.".gene.".$genecounter);

	    if( $ft->has_tag('pseudo') ) {
		$gene->type('pseudo');
	    } else {
		$gene->type('standard');
	    }

	    $trans->id($id.".trans.".$genecounter);
	    $trans->version(1);
	    my $exon = Bio::EnsEMBL::Exon->new();
	    $exon->phase(0);
	    $exon->start($ft->start);
	    $exon->end($ft->end);
	    $exon->strand($ft->strand);
	    $exon->contig_id($self->id);
	    $exon->seqname($self->id);
	    $exon->version(1);
	    $exon->created($time);
	    $exon->modified($time);
	    $exon->id($id.".exon.".$exoncounter++);
	    $trans->add_Exon($exon);

	    my $tranl = Bio::EnsEMBL::Translation->new();
	    $tranl->id($id.".transl.".$genecounter);
	    $tranl->start_exon_id($exon->id);
	    $tranl->end_exon_id($exon->id);
	    $tranl->start(1);
	    $tranl->end($exon->length);
	    $trans->translation($tranl);
	    $tranl->version(1);
	    $genecounter++;
	    push(@genes,$gene);
	} else {
	    # do nothing!
	}
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
   my $length=$self->_get_Seq->length;
   return $length;

}






sub _get_Seq {
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


=head2 embl_offset

 Title   : embl_offset
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub embl_offset{
   my ($self,@args) = @_;

   return 1;
}

=head2 embl_order

 Title   : embl_order
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub embl_order{
   my ($self,@args) = @_;

   return 1;
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
