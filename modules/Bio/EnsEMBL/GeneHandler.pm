
#
# BioPerl module for Bio::EnsEMBL::GeneHandler
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::GeneHandler - Handler for gene objects from ensembl

=head1 SYNOPSIS

    $gh = Bio::EnsEMBL::GeneHandler( -clone => $clone, -gene => $gene);

    # $gh is now a SeqFeature object, with features in the coordinates
    # of the clone

=head1 DESCRIPTION

GeneHandler objects provide the necessary information to 'pretend' that
the stored genes are actually 'real' SeqFeature objects, even though
they aren't.

It provides way to get the 'correct' embl type information out the 
of the objects, potentially with context wrt to the dna sequence.

This is a difficult module to understand as you have to understand how
the AnnSeqIO/EMBL dumper works. It is there where the FtHelper objects
are made.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::GeneHandler;
use vars qw($AUTOLOAD @ISA);
use strict;
use Bio::SeqFeatureI;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;

@ISA = qw(Bio::Root::Object Bio::SeqFeatureI);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;
  my ($clone,$gene) = $self->_rearrange([qw(
					    CLONE
					    GENE )],@args);

  # ok

  $clone || $self->throw("Cannot make a gene handler without a clone");
  $gene  || $self->throw("Cannot make a gene handler without a gene");


  $clone->isa("Bio::EnsEMBL::DB::CloneI") || $self->throw("You haven't given me a valid clone object for this gene!");
  $gene->isa("Bio::EnsEMBL::Gene") || $self->throw("You haven't given me a valid gene object for this gene!");

  
  # this could be a bad way to do this... ;)

  # print STDERR "Dumping in GeneHandler!\n";
  # $gene->_dump(\*STDERR);

  $self->clone($clone);
  $self->gene($gene);
  return $make; # success - we hope!
}

=head2 to_FTHelper

 Title   : to_FTHelper
 Usage   : Provides the necessary transactions for making FTHelper
         : objects for EMBL and GenBank. **Should not be called directly**
         : but only by the embl/genbank feature table dumpers
 Function: 
 Example :
 Returns : 
 Args    :


=cut

sub to_FTHelper{
   my ($self) = @_;
   my (@out);
   my (@exons);

   my $exh = {};
   
   foreach my $trans ( $self->gene()->each_Transcript() ) {
       foreach my $ptrans ( $trans->split_Transcript_to_Partial ) {
	   push(@out,$self->_process_Transcript($ptrans,$exh));
       }
   }


   push(@out,values %{$exh});
   return @out;

}


=head2 clone

 Title   : clone
 Usage   : $obj->clone($newval)
 Function: 
 Returns : value of clone
 Args    : newvalue (optional)


=cut

sub clone{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'clone'} = $value;
    }
    return $obj->{'clone'};

}

=head2 gene

 Title   : gene
 Usage   : $obj->gene($newval)
 Function: 
 Returns : value of gene
 Args    : newvalue (optional)


=cut

sub gene{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'gene'} = $value;
    }
    return $obj->{'gene'};

}


=head2 _process_Transcript

 Title   : _process_Transcript
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _process_Transcript{
   my ($self,$trans,$exon_hash_ref) = @_;

   my $ori;
   my $firstexon;
   my $trans_loc;

   #
   # First exon describes the orientation of the gene
   #
   my @fth;

   ($firstexon) = $trans->each_Exon();
   my $contig = $self->clone()->get_Contig($firstexon->contig_id());
   

   if( $contig->orientation == 1 ) {
       if( $firstexon->strand == 1 ) {
	   $ori = 1;
       } else {
	   $ori = -1;
       }
   } else {
       if( $firstexon->strand == 1 ) {
	   $ori = -1;
       } else {
	   $ori = 1;
       }
   }


   # ok. Now - get into the major loop, and start processing
   # the exons

   foreach my $exon ( $trans->each_Exon() ) {
       # get out the clone
       if( $exon->clone_id() ne $self->clone->id() ) {
	   $self->throw("Cannot currently dump exons across clones");
       }

       $contig = $self->clone()->get_Contig($exon->contig_id());
       $contig->isa("Bio::EnsEMBL::DB::ContigI") || $self->throw("Expecting to get a conting. Instead got a $contig. Not ideal!");

	          
       my ($locstart,$locend,$loc_comp) = $self->_deduce_exon_location($exon,$contig);

       if( ! $exon_hash_ref->{$exon->id()}  ) {
	   # add this exon
	   # make an Exon FTHelper and add them
	   
	   my $ft = new Bio::AnnSeqIO::FTHelper->new();
	   $ft->key("exon");
	   # add other stuff to Exon?
	   $ft->add_field('created',$exon->created());
	   $ft->add_field('modified',$exon->modified());
	   $ft->add_field('exon_id',$exon->id());
	   $ft->add_field('phase',$exon->phase());
	   # $ft->add_field('end_phase',$exon->end_phase());


	   if( $loc_comp == -1 ) {
	       $ft->loc("complement($locstart..$locend)");
	   } else {
	       $ft->loc("$locstart..$locend");
	   }
	   $exon_hash_ref->{$exon->id()} = $ft;
       }

       # add the information to the Transcript object whatever.

       if( $loc_comp == 1 ) {
	   if( ! defined $trans_loc ) {
	       $trans_loc = "<$locstart..$locend";
	   } else {
	       $trans_loc .= ",$locstart..$locend";
	   }
       } else {
	   if( ! defined $trans_loc ) {
	       $trans_loc = "complement(<$locstart..$locend)";
	   } else {
	       $trans_loc .= ",complement($locstart..$locend)";
	   }
       }

       
   }

   # the last *f^%$%ing* location needs to have a '>' (would you
   # believe it). So annoying. This is HORRIBLE

   # puts a > on the last location line
   $trans_loc =~ s/\.\.(\d+)(\)?)$/..$1>$2/;

   my $t_fth = new Bio::AnnSeqIO::FTHelper->new();
   $t_fth->key("CDS");
   my $pseq = $trans->translate();
   $t_fth->add_field('translation',$pseq->str);
   $t_fth->add_field('transcript_id',$trans->id());
   $t_fth->add_field('gene_id',$self->gene->id());
   if( $trans->is_partial() == 1 ) {
       $t_fth->add_field('note',"transcript is a partial transcript");
   }
   
   $t_fth->loc("join($trans_loc)");

   push(@fth,$t_fth);
   return @fth;
}

=head2 _deduce_exon_location

 Title   : _deduce_exon_location
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _deduce_exon_location{
   my ($self,$exon,$contig) = @_;

   my ($locstart,$locend,$loc_comp);

   if( $contig->orientation == 1 ) {
       $loc_comp = $exon->strand;
       if( $loc_comp == 1 ) {
	   $locstart = $contig->offset + $exon->start -1;
	   $locend   = $contig->offset + $exon->end -1;
       } else {
	   $locstart = $contig->offset + $exon->start -1;
	   $locend   = $contig->offset + $exon->end -1;
       }
   } else {
       my $tseq = $contig->seq(); # this is bad news to get out the length
       my $pos  = $contig->offset();

       if( $exon->strand == -1 ) {
	   $locstart = $pos-1 + ($tseq->seq_len() - $exon->end +1);
	   $locend   = $pos-1 + ($tseq->seq_len() - $exon->start +1);
	   $loc_comp = 1;
       } else {
	   $locstart   = $pos-1 + ($tseq->seq_len() - $exon->end +1);
	   $locend     = $pos-1 + ($tseq->seq_len() - $exon->start +1);
	   $loc_comp = -1;
       }
   }
   return ($locstart,$locend,$loc_comp);
}

1;



















