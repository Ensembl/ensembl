
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


=head2 primary_seq

 Title   : seq
 Usage   : $seq = $contig->primary_seq();
 Function: Gets a Bio::PrimarySeqI object out from the contig
 Example :
 Returns : Bio::PrimarySeqI object
 Args    :


=cut

sub pimary_seq {
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

=head1 Decorating methods

These methods do not have to implemented by the derived object.
They are work on top of the interface defined above

=cut


=head2 seq

 Title   : seq
 Usage   : $annseq = $contig->seq();
 Function: Gets a Bio::Seq which can be used as standard a
           Bio::Seq object complete with sequence features (eg genes)
 Example :
 Returns : 
 Args    : Has ref can have a number of attributes to control
           how this call is used.

           $hash->{'strict_EMBL'} = 1; #

           causes EMBL dumping with only EMBL-
           allowed feature qualfiers if true.  Set when
           generating files for submission to the EMBL database.

           $hash->{'seqfeature_filter'} = \&feature_filter_function;

           provides a filter on the sequence features which are attached

=cut

my $global_warn = 0;

sub seq {
    my ($self, $hash_ref) = @_;
    my (@contigs,@genes,$as,$seq);

    if( $global_warn == 1 ) {
	$self->throw("Bad - reentering seq call");
    }

    $global_warn = 1;

    if( defined $hash_ref && ! ref $hash_ref ) {
	$self->throw("Semantics to parameterisation of get_AnnSeq has changed - now pass in a hash");
    }

    if( ! defined $hash_ref ) {
	$hash_ref = {};
    }

    #print STDERR "Starting on the annseq build\n";

    @genes = $self->get_all_Genes();

    #print STDERR "Built genes\n";
    
    $seq = $self->primary_seq();
    
    $as = Bio::Seq->new();
    $as->primary_seq($seq);
    
    # we expect the contigI object to know what id to give ;)
    $as->id($self->id());
    if( $self->can('embl_version') && defined $self->embl_verison ) {
	$as->sv($self->embl_version());
    }
    if( $self->can('htg_phase') && defined $self->htg_phase ) {
	$as->htg_phase($self->htg_phase());
    }

    if( $self->can('seq_date') && defined $self->seq_date ) {
	my $str = POSIX::strftime( "%d-%B-%Y", gmtime($self->seq_date) );
	$as->add_date($str);
    }



    foreach my $gene ( @genes ) {
	print STDERR "got a $gene\n";
        my $gh = new Bio::EnsEMBL::GeneHandler( -gene => $gene,
                                                -strict_embl => $hash_ref->{'strict_EMBL'},
                                                );
        $as->add_SeqFeature($gh);
    }

    #print STDERR "Attached genes\n";

    foreach my $feature ($self->get_all_RepeatFeatures ) {
	$as->add_SeqFeature( $feature );
    }

    return $as;

}


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
