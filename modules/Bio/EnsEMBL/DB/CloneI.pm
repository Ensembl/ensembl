
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
use vars qw($AUTOLOAD @ISA);
use strict;
use Bio::EnsEMBL::GeneHandler;

# Object preamble - inheriets from Bio::Root::Object


@ISA = qw();
# new() is inherited from Bio::Root::Object


=head2 get_Contig

 Title   : get_Contig
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Contig{
   my ($self,@args) = @_;

   $self->warn("Base class has not implemented this yet!");

}


=head2 get_all_Contigs

 Title   : get_all_Contigs
 Usage   : 
 Function: gets all the contigs in this clone
 Example :
 Returns : 
 Args    :


=cut

sub get_all_Contigs{
   my ($self) = @_;

   $self->warn("Base class has not implemented this yet!");

}

=head2 get_all_Genes

 Title   : get_all_Genes
 Usage   :
 Function: gets all the genes that overlap with this clone
 Example :
 Returns : 
 Args    :


=cut

sub get_all_Genes{
   my ($self) = @_;

   $self->warn("Base class has not implemented this yet!");
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

sub seq{
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
	   $self->throw("I am sorry - we have a clone whoses offsets of contigs do not make sense! Contig" . $contig->id() . "at $current_end");
       }

       $seqstr .= 'N' x $nlen;
       
       my $seq = $contig->seq(); # throw exception if it can't do this.
       if( $contig->orientation == -1 ) {
	   $seq = $seq->revcom();
       }

       $seqstr .= $seq->str();

       $current_end = $contig->offset + $seq->seq_len;
   }

   $out = Bio::Seq->new( -id => $self->id() , -seq => $seqstr, -type => 'Dna');

   return $out;
}

=head2 get_AnnSeq

 Title   : get_AnnSeq
 Usage   : $annseq = $clone->get_AnnSeq()
 Function: Gets a Bio::AnnSeq which can be used as standard
 Example :
 Returns : 
 Args    :


=cut

sub get_AnnSeq{
   my ($self) = @_;

   my (@contigs,@genes,$as,$seq);

   @contigs = $self->get_all_Contigs();
   @genes   = $self->get_all_Genes();

   $seq = $self->seq();
   
   $as = Bio::AnnSeq->new();
   $as->seq($seq);
   foreach my $gene ( @genes ) {
       print STDERR "Adding gene $gene\n";
       my $gh = new Bio::EnsEMBL::GeneHandler( -clone => $self, -gene => $gene );
       $as->add_SeqFeature($gh);
   }

   return $as;
}

1;

