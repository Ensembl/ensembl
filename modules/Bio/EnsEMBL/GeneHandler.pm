
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

   my %exh;
   
   foreach my $trans ( $self->gene()->each_Transcript() ) {
       my $trans_strand = undef;
       my $trans_loc;
       # wrap this in an eval to catch database exceptions.
       eval {
	   foreach my $exon ( $trans->each_Exon() ) {
	       # get out the clone
	       if( $exon->clone_id() eq $self->clone->id() ) {

		   # in the future we probably need to factor out
		   # this inner loop as we will need to reuse it

		   # find the offset position of the contig in the clone
		   my $contig = $self->clone()->get_Contig($exon->contig_id());
		   my $pos = $contig->offset();
		   my $ori = $contig->orientation();
		   
                   # make an Exon FTHelper and add them

		   my $ft = new Bio::AnnSeqIO::FTHelper->new();

		   
		   
		   my $locstart;
		   my $locend;
		   my $exst;
		   
		   if( $ori == 1 ) {
		       $locstart = $exon->start+$pos-1;
		       $locend   = $exon->end+$pos-1;
		       $ft->loc("$locstart..$locend");
		       if( $exon->strand == -1 ) {
			   $exst = -1;
			   $ft->loc("complement(". $ft->loc . ")");
		       } else {
			   $exst = 1;
		       }
		   } else {
		       $locstart = $pos-1 + ($contig->length() - $exon->end +1);
		       $locend   = $pos-1 + ($contig->length() - $exon->start +1);
		       $ft->loc("$locstart..$locend"); 
		       if( $exon->strand == 1 ) {
			   $exst = -1;
			   $ft->loc("complement(". $ft->loc . ")");
		       } else {
			   $exst = 1;
		       }
		   }
		   $ft->key("Exon");

		   # add other stuff to Exon?
		   $ft->add_field('created',$exon->created());
		   $ft->add_field('modified',$exon->modified());
		   $ft->add_field('exon_id',$exon->id());
		   $ft->add_field('phase',$exon->phase());
		   $ft->add_field('end_phase',$exon->end_phase());

		   # FIXME: we are redoing too much information above
                   # add this to a hash so that we can get out the 
                   # unique set of exons for this gene. We can't stop
                   # processing this exon, as we need the information for transcript.

		   $exh{$exon->id()} = $ft;
		   

		   # ok - now handle the location line in the transcript object
		   if( defined $trans_loc  ) {
		       $trans_loc .= ",";
		       $trans_loc .= "$locstart..$locend";
		   } else {
		       $trans_loc = "";
		       $trans_loc = "$locstart..$locend";
		   }


		   if( !defined $trans_strand ) {
		       $trans_strand = $exst;
		   } else {
		       if( $exst != $trans_strand ) {
			   $self->warn("Oh no - transcript with orientation different from implied gene order");
		       }
		   }

		   
	       } else {
		   $self->throw("Have not delt with exons on other clones yet! Self is ". $self->clone->id(). " exon is " . $exon->clone_id());
	       }
	   }

	   my $fth = Bio::AnnSeqIO::FTHelper->new();
	   $fth->key("CDS");
	   
	   $trans_loc = "join($trans_loc)";
	   if( $trans_strand == -1 ) {
	       $fth->loc("complement($trans_loc)");
	   } else {
	       $fth->loc($trans_loc);
	   }
	   my $pseq = $trans->translate();
	   $fth->add_field('translate',$pseq->str);
	   $fth->add_field('transcript_id',$trans->id());
	   $fth->add_field('gene_id',$self->gene->id());
	   
	   push(@out,$fth);
       };
       if( $@ ) {
	   $self->throw("Unable to process transcript due to\n$@\n");
       }
   }
   push(@out,values %exh);
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

