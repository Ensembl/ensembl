
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


    # contigs can be made in a number of different ways
    $contig = $obj->get_Contig($contigid);
 
    # contigs objects have an extend method. This gives back
    # a new contig object which is longer to the 5' and 3' end
    # If it runs out of sequence, it truncates silently.
    $virtual_contig = $contig->extend(1000,1000);

    # contigs have special feature extraction functions
    @repeats = $contig->get_all_RepeatFeatures();
    @sim     = $contig->get_all_SimilarityFeatures();

    # you can get genes attached to this contig. This does not
    # mean that all the gene is on this contig, just one exon
    @genes   = $contig->get_all_Genes();

    # ContigI is-a Bio::SeqI which is-a PrimarySeqI. This means
    # that the normal bioperl functions work ok. For example:

    $string = $contig->seq(); # the entire sequence
    $string = $contig->subseq(100,120);  # a sub sequence

    $seqout = Bio::SeqIO->new( '-format' => 'embl', -fh => \*STDOUT );
    $seqout->write_seq($contig);


=head1 DESCRIPTION

The contig interface defines a single continuous piece of DNA with both
features and genes on it. It is-a Bio::SeqI interface, meaning that it
can be used in any function call which takes bioperl Bio::SeqI objects.

It has additional methods, in particular the ability to only get a 
subset of features out and genes.

The contig interface just defines a number of functions which have to provided 
by implementations. Two good implementations are the RawContig implementation
found in Bio::EnsEMBL::DBSQL::RawContig and the generic VirtualContig interface
in Bio::EnsEMBL::DB::VirtualContig

=head1 CONTACT

Ewan Birney, <birney@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DB::ContigI;

use vars ('@ISA');
use strict;
use Bio::EnsEMBL::VirtualGene;
use Bio::SeqI;
use Bio::Root::RootI

@ISA = qw( Bio::SeqI Bio::Root::RootI );

=head2 primary_seq

 Title   : seq
 Usage   : $seq = $contig->primary_seq();
 Function: Gets a Bio::PrimarySeqI object out from the contig
 Example :
 Returns : Bio::PrimarySeqI object
 Args    :


=cut

sub primary_seq {
   my ($self) = @_;
   $self->throw("Object did not provide the primary_seq method on a contig interface");
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

   $self->throw("Object did not provide the get_all_SimilarityFeatures method on Contig interface!");

}

=head2 get_all_RepeatFeatures

 Title   : get_all_RepeatFeatures
 Usage   : foreach my $sf ( $contig->get_all_RepeatFeatures ) 
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_RepeatFeatures{
   my ($self) = @_;

   $self->throw("Object did not provide the get_all_RepeatFeatures method on Contig interface!");

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

=head2 extend

 Title   : extend
 Usage   : $newcontig = $contig->extend(1000,-1000)
 Function: Makes a new contig shifted along by the base pairs to the
           5' and the 3'. 
 Example :
 Returns : A ContigI implementing object
 Args    :


=cut

sub extend{
   my ($self,@args) = @_;

   $self->throw("Object did not provide the extend method on Contig interface!");
}

=head2 dbobj

 Title   : dbobj
 Usage   : $obj = $contig->dbobj
 Function: returns a Bio::EnsEMBL::DB::ObjI implementing function
 Example :
 Returns : 
 Args    :


=cut

sub dbobj{
   my ($self,@args) = @_;

   $self->throw("Object did not provide the dbobj method on the Contig interface");
}


=head2 SeqI implementing methods

As ContigI is-a SeqI, we need to implement some sequence
feature methods. It is in these calls where the "magic" happens
by calling VirtualGene for genes to map genes to contigs.

You do not need to implement this methods, but you can if you 
wish to control their behaviour

=cut

=head2 top_SeqFeatures

 Title   : top_SeqFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub top_SeqFeatures{
   my ($self,@args) = @_;
   my (@f);
   push(@f,$self->get_all_SeqFeatures());
   foreach my $gene ( $self->get_all_Genes()) {
       print STDERR "Got a $gene\n";
       my $vg = Bio::EnsEMBL::VirtualGene->new(-gene => $gene,-contig => $self);
       push(@f,$vg);
   }

   return @f;
}



=head2 all_SeqFeatures

 Title   : all_SeqFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub all_SeqFeatures {
   my ($self) = @_;
   my (@array);
   foreach my $feat ( $self->top_SeqFeatures() ){
       push(@array,$feat);
       &_retrieve_subSeqFeature(\@array,$feat);
   }

   return @array;
}


sub _retrieve_subSeqFeature {
    my ($arrayref,$feat) = @_;

    foreach my $sub ( $feat->sub_SeqFeature() ) {
	push(@$arrayref,$sub);
	&_retrieve_subSeqFeature($arrayref,$sub);
    }

}

=head2 annotation

 Title   : annotation
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub annotation{
   my ($self,@args) = @_;

   if( $self->can('get_annotation_hook') ) {
       return $self->get_annotation_hook();
   }
   return ();
}


=head2 PrimarySeqI implementing methods

As Bio::SeqI is-a PrimarySeqI, we need to implement these methods.
They can all be delegated to PrimarySeq. You do not need to implement
these methods

=cut

=head2 seq

 Title   : seq
 Usage   : $string = $contig->seq();
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub seq{
   my ($self,@args) = @_;

   return $self->primary_seq->seq();
}

=head2 subseq

 Title   : subseq
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub subseq{
   my ($self,$start,$end) = @_;

   return $self->primary_seq->subseq($start,$end);

}

=head2 display_id

 Title   : display_id
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub display_id{
   my ($self,@args) = @_;

   return $self->id();
}

=head2 primary_id

 Title   : primary_id
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub primary_id{
   my ($self,@args) = @_;

   return "$self";
}

=head2 accession_number

 Title   : accession_number
 Usage   : $obj->accession_number($newval)
 Function: 
 Returns : value of accession_number
 Args    : newvalue (optional)


=cut

sub accession_number{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'accession_number'} = $value;
    }
    return $obj->{'accession_number'};

}

=head2 desc

 Title   : desc
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub desc{
   my ($self,@args) = @_;

   return "Ensembl Contig";
}


=head1 Decorating methods

These methods do not have to implemented by the derived object.
They are work on top of the interface defined above

=cut



sub get_AnnSeq {
    my $self = shift;
    $self->throw("You should use seq function on the ContigI interface");
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


=head1 Cruft

Not clear if this method belongs here....

=cut

#
# Not sure where to put this?
#
 
=head2 find_supporting_evidence

 Title   : find_supporting_evidence
 Usage   : $obj->find_supporting_evidence($exon);
 Function: Looks through all the similarity features and
           stores as supporting evidence any feature
           that overlaps with an exon.  I know it is
           a little crude but it\'s a start/
 Example : 
 Returns : Nothing
 Args    : Bio::EnsEMBL::Exon


=cut


sub find_supporting_evidence {
    my ($self,$exon) = @_;

    my @features = $self->get_all_SimilarityFeatures;

    foreach my $f (@features) {
	if ($f->overlaps($exon)) {
	    $exon->add_Supporting_Feature($f);
	}
    }
}
    

1;
