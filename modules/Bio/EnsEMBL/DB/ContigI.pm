
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


package Bio::EnsEMBL::DB::ContigI;

use vars qw(@ISA);
use strict;
use Bio::AnnSeq;

=head2 seq

 Title   : seq
 Usage   : $seq = $contig->seq();
 Function: Gets a Bio::Seq object out from the contig
 Example :
 Returns : Bio::Seq object
 Args    :


=cut

sub seq{
   my ($self) = @_;
   $self->throw("Object did not provide the seq method on a contig interface");
}

=head2 get_all_SeqFeatures

 Title   : get_all_SeqFeatures
 Usage   : foreach my $sf ( $contig->get_all_SeqFeatures ) 
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_SeqFeatures{
   my ($self) = @_;

   $self->throw("Object did not provide the get_all_SeqFeatures method on Contig interface!");

}


=head2 get_all_Genes

 Title   : get_all_Genes
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_Genes{
   my ($self) = @_;

   $self->throw("Object did not provide the get_all_Genes method on Contig interface!");

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

   $self->throw("Object did not provide the offset method on Contig interface!");

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

   $self->throw("Object did not provide the orientation method on Contig interface!");

}

=head2 annseq

 Title   : annseq
 Usage   : $annseq = $contig->annseq();
 Function: Gets Bio::AnnSeq object, for example, 
           ready for EMBL dumping
 Returns : Bio::AnnSeq object
 Args    :


=cut

sub annseq{
   my ($self) = @_;

   my $seq = $self->seq();
   my $annseq = Bio::AnnSeq->new();
   
   foreach my $sf ( $self->get_all_SeqFeatures() ) {
       $annseq->add_SeqFeature($sf);
   }
   $annseq->seq($seq);

   return $annseq;

}


=head2 write_acedb

 Title   : write_acedb
 Usage   : $contig->write_acedb(\*FILEHANDLE);
           $contig->write_acedb(\*FILEHANDLE,$ace_seq_name);
 Function: Dumps exon, transcript and gene objects in acedb format
 Returns : 
 Args    :

=cut

sub write_acedb{
    my ($self,$fh,$seqname) = @_;

    my $contig_id=$self->id();

    $seqname ||= $contig_id;
    
    foreach my $gene ($self->get_all_Genes()){
	my $gene_id=$gene->id;
	TRANSCRIPT :
	foreach my $trans ( $gene->each_Transcript ) {
	    my $trans_id=$trans->id;
	    
	    # check this transcript has exons on this contig
	    foreach my $exon ( $trans->each_Exon ) {
		if( $exon->contig_id ne $contig_id ) {
		    $self->warn("Could not ace dump transcript " . $trans->id . "as exons across contigs");
		    next TRANSCRIPT;
		}
	    }
	    
	    # exons are in order.

	    my @exons = $trans->each_Exon;

	    my $tstrand = $exons[0]->strand;
	    my ($tstart,$tend);
	    if( $tstrand == 1 ) {
		$tstart = $exons[0]->start;
		$tend   = $exons[$#exons]->end;
	    } else {
		$tstart = $exons[0]->end;
		$tend   = $exons[$#exons]->start;
	    }

	    # print starting stuff...

	    print $fh "Sequence $seqname\n";
	    print $fh "subsequence $gene_id.$trans_id.EnsEMBL $tstart $tend\n\n";
		
	    # acedb has coordinates relative to transcripts.
	    
	    print $fh "Sequence $gene_id.$trans_id.EnsEMBL\nCDS\nStart_not_found\nEnd_not_found\n";
	    
	    foreach my $exon ( $trans->each_Exon ) {
		if( $tstrand == 1 ) {
		    print $fh "source_Exons ", ($exon->start - $tstart + 1)," ",($exon->end - $tstart +1), "\n";
		} else {
		    print $fh "source_Exons ", ($tstart - $exon->end +1 ), " ",($tstart - $exon->start+1),"\n";
		}
	    }

	    print $fh "\n\n";
	}
    }

}




=head2 as_seqfeatures

 Title   : as_seqfeatures
 Usage   : @seqfeatures = $contig->as_seqfeatures();
           foreach $sf ( @seqfeatures ) { 
	       print $sf->gff_string(), "\n";
           }
 Function: Makes ensembl exons as an array of seqfeature::generic
           objects that can be dumped with the correct additional tags
           about transcripts/genes etc added to them
 Returns : An array of SeqFeature::Generic objects
 Args    :

=cut

sub as_seqfeatures {
    my ($self) = @_;
    my $contig_id=$self->id();
    my @sf;

    foreach my $gene ($self->get_all_Genes()){
	my $gene_id=$gene->id;
	foreach my $trans ( $gene->each_Transcript ) {
	    my $transcript_id=$trans->id;
	    foreach my $exon ( $trans->each_Exon ) {
		my $sf= Bio::SeqFeature::Generic->new();
		$sf->seqname($contig_id);
		$sf->source_tag('ensembl');
		$sf->primary_tag('exon');
		$sf->start($exon->start);
		$sf->end($exon->end);
		$sf->strand($exon->strand);
		#$sf->frame($exon->frame);
		$sf->add_tag_value('ensembl_exon_id',$exon->id);
		$sf->add_tag_value('ensembl_transcript_id',$transcript_id);
		$sf->add_tag_value('ensembl_gene_id',$gene_id);
		push(@sf,$sf);
	    }
	}

    }
    return @sf;
}


1;
