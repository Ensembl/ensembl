
#
# BioPerl module for Transcript
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Transcript - gene transcript object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Contains details of coordinates of all exons that make
up a gene transcript.

Creation:
   
     my $tran = new Bio::EnsEMBL::Transcript();
     my $tran = new Bio::EnsEMBL::Transcript(@exons);

Manipulation:

     my @exons = $tran->each_Exon         # Returns an array of Exon objects
     my $pep   = $tran->translate()       # Returns the peptide translation of the exons as a Bio::Seq
     
     $tran->sort()                        # Sorts exons into order (forward for + strand, reverse fo - strand)

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::Transcript;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object

use Bio::Root::Object;


@ISA = qw(Bio::Root::Object);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  $self->{'_trans_exon_array'} = [];
  my $make = $self->SUPER::_initialize;

  # set stuff in self from @args
  foreach my $a (@args) {
    $self->add_Exon($a);
  }

  return $self; # success - we hope!
}

=head2 id

 Title   : id
 Usage   : $obj->id($newval)
 Function: 
 Returns : value of id
 Args    : newvalue (optional)


=cut

sub id{
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'id'} = $value;
    }
    return $self->{'id'};

}

=head2 add_Exon

 Title   : add_Exon
 Usage   : $trans->add_Exon($exon)
 Returns : Nothing
 Args    :


=cut

sub add_Exon{
   my ($self,$exon) = @_;

   #yup - we are going to be picky here...

   if( ! $exon->isa("Bio::EnsEMBL::Exon") ) {
       $self->throw("$exon is not a Bio::EnsEMBL::Exon!");
   }

   # at the moment, use the SeqFeature sub hash. But in the future,
   # possibly do something better?

   push(@{$self->{'_trans_exon_array'}},$exon);
   
}

=head2 each_Exon

 Title   : each_Exon
 Usage   : foreach $exon ( $trans->each_Exon)
 Function: Returns an array of exons in the transcript
 Example : my @exons = $tr->each_Exon
 Returns : An array of exon objects
 Args    : none


=cut

sub each_Exon{
   my ($self) = @_;
   my @sub;
   my @ret;

   return @{$self->{'_trans_exon_array'}};
}


=head2 flush_Exon

 Title   : flush_Exon
 Usage   : Removes all Exons from the array.
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub flush_Exon{
   my ($self,@args) = @_;

   $self->{'_trans_exon_array'} = [];
}

=head2 translate

 Title   : translate
 Usage   : $pep = $feat->translate()
 Function: Returns the peptide translation of the gene - in the correct phase
 Returns : Bio::Seq
 Args    : none

=cut

sub translate {
  my ($self) = @_;
  my @exon = $self->each_Exon;
  
  my $phase = $exon[0]->phase;
  my $mrna  = $self->dna_seq->seq();

  # Hackeroover as we can't ranslate in any frame apart from 0
  $mrna = substr($mrna,$phase);

  my $seq = new Bio::Seq(-seq => $mrna);
  my $trans = $seq->translate();

  return $trans;
}

=head2 dna_seq

 Title   : dna_Seq
 Usage   : $dna = $feat->dna_seq
 Function: Returns the dna sequence of the gene
 Returns : Bio::Seq
 Args    : none

=cut

sub dna_seq {
  my ($self) = @_;

  my $mrna = "";

  foreach my $exon ($self->each_Exon) {
    $mrna .= $exon->dna_seq()->seq();
  }

  my $seq = new Bio::Seq(-seq => $mrna);

  return $seq;
}

=head2 sort

 Title   : sort
 Usage   : $feat->sort()
 Function: Sorts the exon features by start coordinate
           Sorts forward for forward strand and reverse for reverse strand
 Returns : none
 Args    : none

=cut

sub sort {
  my $self = shift;

  # Fetch all the features
  my @exons = $self->each_Exon();

  # Empty the feature table
  $self->flush_Exon();

  # Now sort the exons and put back in the feature table
  my $strand = $exons[0]->strand;

  if ($strand == 1) {
    @exons = sort { $a->start <=> $b->start } @exons;
  } elsif ($strand == -1) {
    @exons = sort { $b->start <=> $a->start } @exons;
  }

  foreach my $e (@exons) {
    $self->add_Exon($e);
  }
}


1;
