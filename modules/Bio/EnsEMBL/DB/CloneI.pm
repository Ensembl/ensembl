
#
# BioPerl module for Bio::EnsEMBL::DB::CloneI
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DB::CloneI - Abstract Interface of a clone object

=head1 SYNOPSIS

    # get a clone object somehow

    @contigs = $clone->get_all_Contigs();

    @genes   = $clone->get_all_Genes();

    # dumping EMBL format

    $ostream = Bio::AnnSeqIO->new( -format => 'EMBL' , -fh => \*STDOUT );
    $annseq  = $clone->get_AnnSeq();
    $ostream->write_annseq($annseq);

    
=head1 DESCRIPTION

Ewan Birney

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DB::CloneI;
use vars qw($AUTOLOAD @ISA $CONTIG_SPACING);
use strict;
use Bio::EnsEMBL::GeneHandler;
use Bio::EnsEMBL::AnnSeq;
use POSIX;

# Object preamble - inheriets from Bio::Root::Object

$CONTIG_SPACING = 800;

@ISA = qw();
# new() is inherited from Bio::Root::Object


=head2 id

 Title   : id
 Usage   : this is the primary id for ensembl. General the accession number from embl.
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub id {
   my ($self,@args) = @_;

   $self->warn("Base class has not implemented this yet!");

}

=head2 embl_id

 Title   : embl_id
 Usage   : this is the embl_id for this clone, to generate nice looking files
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub embl_id {
   my ($self,@args) = @_;

   $self->warn("Base class has not implemented embl_id yet!");

}


=head2 sv

 Title   : sv
 Function: returns the version number (not the acc.version, just verision).
 Example :
 Returns : 
 Args    :


=cut

sub sv {
   my ($self,@args) = @_;

   $self->warn("Base class has not implemented this yet!");

}


=head2 embl_version

 Title   : embl_version
 Usage   : $clone->embl_version()
 Function: Gives the value of the EMBL version, i.e. the data version
 Example : $clone->embl_version()
 Returns : version number
 Args    : none


=cut

sub embl_version {
    my ($self,@args) = @_;

    $self->warn("Base class has not implemented embl_version yet!");

}


=head2 seq_date

 Title   : seq_date
 Usage   : $clone->seq_date()
 Function: loops over all $contig->seq_date, throws a warning if they are different and 
           returns the first unix time value of the dna created datetime field, which indicates
           the original time of the dna sequence data
 Example : $clone->seq_date()
 Returns : unix time
 Args    : none


=cut

sub seq_date {
    my ($self,@args) = @_;

    $self->warn("Base class has not implemented seq_date yet!");
}


=head2 version

 Title   : version
 Function: Schema translation
 Example :
 Returns : 
 Args    :


=cut

sub version {
   my ($self,@args) = @_;

   $self->warn("Called version without implementation. Probably an old object. Called sv instead");

   return $self->sv();
}


=head2 htg_phase

 Title   : htg_phase
 Usage   : this is the phase being 1,2,3,4 (4 being finished).
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub htg_phase {
   my ($self,@args) = @_;

   $self->warn("Base class has not implemented htg_phase yet!");

}


=head2 created

 Title   : created
 Usage   : $clone->created()
 Function: Gives the unix time value of the created datetime field, which indicates
           the first time this clone was put in ensembl
 Example : $clone->created()
 Returns : unix time
 Args    : none


=cut

sub created {
    my ($self,@args) = @_;
    
    $self->warn("Base class has not implemented created  yet!");
}


=head2 modified

 Title   : modified
 Usage   : $clone->modified()
 Function: Gives the unix time value of the modified datetime field, which indicates
           the last time this clone was modified in ensembl
 Example : $clone->modified()
 Returns : unix time
 Args    : none


=cut

sub modified{
    my ($self,@args) = @_;

    $self->warn("Base class has not implemented modified yet!");
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
   my ($self,@args) = @_;

   $self->warn("Base class has not implemented get_Contig yet!");

}


=head2 get_all_Contigs

 Title   : get_all_Contigs
 Usage   : 
 Function: gets all the contigs in this clone
 Example :
 Returns : 
 Args    :


=cut

sub get_all_Contigs {
   my ($self) = @_;

   $self->warn("Base class has not implemented get_all_Contigs yet!");

}


=head2 get_all_Genes

 Title   : get_all_Genes
 Usage   :
 Function: gets all the genes that overlap with this clone
 Example :
 Returns : 
 Args    :


=cut

sub get_all_Genes {
   my ($self) = @_;

   $self->warn("Base class has not implemented get_all_Genes yet!");
}


=head2 seq

 Title   : seq
 Usage   : $seq = $clone->seq()
 Function: Provides a Bio::Seq object which represents the
           the clone, potentially with N's inserted. 
 Example :
 Returns : A Bio::Seq object
 Args    : 


=cut

sub seq {
   my ($self,$spacer) = @_;
   my $out;
   my $seqstr = "";
   my $current_end;
   
   
   my @contigs = $self->get_all_Contigs();
   # get paranoid about contigs with their orientation!
   
   @contigs = sort { $a->offset <=> $b->offset } @contigs;
   $current_end = 1;
   foreach my $contig ( $self->get_all_Contigs ) {
       my $nlen = $contig->offset - $current_end;
       if( $nlen < 0 ) {
	   $self->throw("ERROR: offsets of contigs do not make sense! Contig " . 
			$contig->id() . " at $current_end");
       }

       $seqstr .= 'N' x $nlen;
       
       my $seq = $contig->seq(); # throw exception if it can't do this.

       if( ! $seq->isa('Bio::Seq') ) {
	   $self->throw("Got a $seq not a Bio::Seq!");
       }

       #if( $seq->type ne 'Dna' ) {
#	   $self->warn("For contig " . $contig->id . "sequence type is ". $seq->type . " not Dna");
 #      }

       if( $contig->orientation == -1 ) {
	   $seq = $seq->revcom();
       }

       $seqstr .= $seq->str();

       $current_end = $contig->offset + $seq->seq_len;
   }

   $out = Bio::Seq->new( '-id' => $self->id() , -seq => $seqstr, -type => 'Dna');

   return $out;
}


=head2 get_AnnSeq

 Title   : get_AnnSeq
 Usage   : $annseq = $clone->get_AnnSeq($referece_to_hash)
 Function: Gets a Bio::AnnSeq which can be used as standard
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

sub get_AnnSeq {
    my ($self, $hash_ref) = @_;
    my (@contigs,@genes,$as,$seq);

    if( defined $hash_ref && ! ref $hash_ref ) {
	$self->throw("Semantics to parameterisation of get_AnnSeq has changed - now pass in a hash");
    }

    if( ! defined $hash_ref ) {
	$hash_ref = {};
    }

    #print STDERR "Starting on the annseq build\n";

    @genes = $self->get_all_Genes();

    #print STDERR "Built genes\n";
    
    $seq = $self->seq();
    
    $as = Bio::EnsEMBL::AnnSeq->new();
    
    $as->embl_id($self->embl_id());
    $as->sv($self->embl_version());
    $as->htg_phase($self->htg_phase());

    $as->seq($seq);

    my $str = POSIX::strftime( "%d-%b-%Y", gmtime($self->seq_date) );
    $as->add_date($str);



    foreach my $gene ( @genes ) {
        my $gh = new Bio::EnsEMBL::GeneHandler( -clone => $self,
                                                -gene => $gene,
                                                -strict_embl => $hash_ref->{'strict_EMBL'},
                                                );
        $as->add_SeqFeature($gh);
    }

    #print STDERR "Attached genes\n";

    # Add features to annseq object
    foreach my $contig ($self->get_all_Contigs) {
        print STDERR "Getting features for contig" .$contig->id."\n";

        # Coordinates retrieved are in Clone coordinate space
        # from the get_all_clone_SeqFeatures method

	#
        foreach my $feature ($contig->get_clone_RepeatFeatures ) {

	    # apply filter if there
	    if( $hash_ref->{'seqfeature_filter'} ) {
		if( &{$hash_ref->{'seqfeature_filter'}}($feature) != 1 ) {
		    next;
		}
	    }
            $as->add_SeqFeature( $feature );
        }
    }
    print STDERR "Built AnnSeq\n";
    return $as;
}

1;
