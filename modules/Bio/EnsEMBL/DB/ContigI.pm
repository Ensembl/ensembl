
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

1;
