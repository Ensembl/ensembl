
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

use strict;
use Bio::AnnSeq;


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
   $self->throw("Class [$self] has not implemented the created method");
}

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

=head2 id

 Title   : id
 Usage   : $obj->id($newval)
 Function: 
 Example : 
 Returns : value of id
 Args    : newvalue (optional)


=cut

sub id{
    my ($self) = @_;
    $self->throw("Object did not provide the id method on a contig interface");
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

=head2 get_all_clone_SeqFeatures

 Title   : get_all_clone_SeqFeatures
 Usage   : foreach $feat ( $contig->get_all_clone_SeqFeatures )
 Function: returns sequence features but in the clone coordinate space.
 Example :
 Returns : an array of SeqFeatures
 Args    : None

 This method is common to all Contig objects, whatever the implementation.
Implementation objects do not need to write this method

=cut

sub get_all_clone_SeqFeatures{
   my ($self) = @_;
   my @out;

   foreach my $sf ( $self->get_all_SeqFeatures ) {
       my ($start,$end,$strand) = $self->_convert_coords_contig_clone($sf->start,$sf->end,$sf->strand);
       $sf->start($start);
       $sf->end($end);
       $sf->strand($strand);
       push(@out,$sf);
   }

   return @out;
}


sub _convert_coords_contig_clone {
    my $self = shift;
    my $start = shift;
    my $end = shift;
    my $strand = shift;

    my ($out_start,$out_end,$out_strand);

    if( !defined $strand ) {
	$self->throw("_convert_coords_contig_clone(start,end,strand)");
    }

    my $offset = $self->offset;

    if( $self->orientation == 1 ) {
       $out_strand = $strand;
       if( $out_strand == 1 ) {
	   $out_start = $offset + $start -1;
	   $out_end   = $offset + $end -1;
       } else {
	   $out_start = $offset + $start -1;
	   $out_end   = $offset + $end -1;
       }
   } else {
       my $length = $self->length(); 
       if( $strand == -1 ) {
	   $out_start = $offset-1 + ($length - $end +1);
	   $out_end   = $offset-1 + ($length - $start +1);
	   $out_strand = 1;
       } else {
	   $out_start   = $offset-1 + ($length - $end +1);
	   $out_end     = $offset-1 + ($length - $start +1);
	   $out_strand = -1;
       }
   }

   return ($out_start,$out_end,$out_strand);
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

    $self->throw("Object did not provide the seq_date method on Contig interface!");

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

=head2 order

 Title   : order
 Usage   : $obj->order($newval)
 Function: 
 Returns : value of order
 Args    : newvalue (optional)


=cut

sub order{
    my ($self,@args) = @_;

    $self->throw("Object did not provide the order method on Contig interface!");
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

   $self->throw("Object did not provide the length method on Contig interface!");

}

=head2 annseq

 Title   : annseq
 Usage   : $annseq = $contig->annseq();
 Function: Gets Bio::AnnSeq object, for example, 
           ready for EMBL dumping
 Returns : Bio::AnnSeq object
 Args    :


=cut

sub annseq {
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

sub write_acedb {
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
	    
	    print $fh "Sequence $gene_id.$trans_id.EnsEMBL\nCDS\nStart_not_found\nEnd_not_found\nMethod EnsEMBL\n";
	    
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

    # build objects for each exon in each gene
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
		$sf->add_tag_value('contig_id',$contig_id);
		push(@sf,$sf);
	    }
	}
    }

    # add objects for each feature on contig
    push(@sf,$self->get_all_SeqFeatures);

    return @sf;
}


1;
