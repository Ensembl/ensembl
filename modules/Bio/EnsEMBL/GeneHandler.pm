
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
       push(@out,$self->_process_Transcript($trans,$exh));
   }
    #   foreach my $ptrans ( $trans->split_Transcript_to_Partial ) {
    #	   push(@out,$self->_process_Transcript($ptrans,$exh));
    #   }



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

   # we want to know whether the last exon is forward or backward first.
   my $loc_comp; 
   my $prev;
   print STDERR "Looking at this transcript ",$trans->id(),"\n";

   foreach my $exon ( $trans->each_Exon() ) {
       # get out the clone

       print STDERR "Looking at exon ",$exon->id,"\n";

       if( $exon->clone_id() ne $self->clone->id() ) {
	   $self->throw("Cannot currently dump exons across clones");
       }

       $contig = $self->clone()->get_Contig($exon->contig_id());
       $contig->isa("Bio::EnsEMBL::DB::ContigI") || $self->throw("Expecting to get a conting. Instead got a $contig. Not ideal!");

       my ($locstart,$locend);
       ($locstart,$locend,$loc_comp) = $self->_deduce_exon_location($exon,$contig);
       

       if( ! $exon_hash_ref->{$exon->id()}  ) {
	   # add this exon
	   # make an Exon FTHelper and add them
	   
	   my $ft = new Bio::AnnSeqIO::FTHelper->new();
	   $ft->key("exon");
	   # add other stuff to Exon?
	   $ft->add_field('created',$exon->created());
	   $ft->add_field('modified',$exon->modified());
	   $ft->add_field('exon_id',$exon->id());
	   $ft->add_field('start_phase',$exon->phase());
	   $ft->add_field('end_phase',$exon->end_phase());


	   if( $loc_comp == -1 ) {
	       $ft->loc("complement($locstart..$locend)");
	   } else {
	       $ft->loc("$locstart..$locend");
	   }
	   $exon_hash_ref->{$exon->id()} = $ft;
       }

       # check to see whether we are jumping contigs or not

       if( $prev && $prev->contig_id ne $exon->contig_id ) {
	   my $loc_exon = $self->_generate_missing_exon($self->clone->get_Contig($prev->contig_id),
						 $self->clone->get_Contig($exon->contig_id),
						 (3-$prev->end_phase) + (3-((3-$exon->phase)%3)),$exon->strand); 
	   # generate missing exon now, in CDS line

	   print STDERR "Adding $loc_exon\n";

	   $trans_loc .= ",$loc_exon";
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
	       $trans_loc = "complement($locstart..>$locend)";
	   } else {
	       $trans_loc .= ",complement($locstart..$locend)";
	   }
       }

       $prev = $exon;
   } # end of loop over all exons

   # the last *f^%$%ing* location needs to have a '>' (would you
   # believe it). So annoying. This is ..HORRIBLE..


   # FIXME:
   # better solution to deal with the last exon separately

   # puts a > on the last location line

   if( $loc_comp == 1 ) {
       $trans_loc =~ s/(\d+)\.\.(\d+)(\)?)$/$1..>$2$3/;
   } else {
       $trans_loc =~ s/(\d+)\.\.(>?)(\d+)(\)?)$/<$1..$2$3$4/;
   }

   # use first exon to find appropiate start point


   my $t_fth = new Bio::AnnSeqIO::FTHelper->new();
   $t_fth->key("CDS");
   my $pseq = $trans->translate();
   $t_fth->add_field('translation',$pseq->str);
   $t_fth->add_field('transcript_id',$trans->id());
   $t_fth->add_field('gene_id',$self->gene->id());

   # map phase1 -> 2, phase 2 -> 1, phase 0 -> 0 and then add 1. Easy huh?
   $t_fth->add_field('codon_start',((3 -$firstexon->phase) %3)+1);

   if( $trans->is_partial() == 1 ) {
       $t_fth->add_field('note',"transcript is a partial transcript");
   }

   # hacky way of figuring out whether we need to "join" or not
   if( $trans_loc =~ /,/ ) {
       $t_fth->loc("join($trans_loc)");
   } else {
       $t_fth->loc($trans_loc);
   }


   push(@fth,$t_fth);
   return @fth;
}

=head2 _generate_missing_exon

 Title   : _generate_missing_exon
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _generate_missing_exon{
   my ($self,$contiga,$contigb,$number,$strand) = @_;

   if( ! $strand || !$contiga->isa('Bio::EnsEMBL::DB::ContigI') || !$contigb->isa('Bio::EnsEMBL::DB::ContigI') ){
       $self->throw("bad arguments, $contiga $contigb");
   }

   # yuk yuk yuk.
  
   # issue a warning if contiga is not next door to contigb.

   if( abs($contiga->order - $contigb->order) != 1 ) {
       $self->warn("exons spanning more than one gap - ie a contig.");
   }

   # take contigb's offset - take mid point of joining segment to contigb

   my $start;
   my $end;
   my $isc;

   # figure out whether conitga is before or after contigb
   my $corder;
   if( $contiga->order < $contigb ) {
       $corder = 1;
   } else {
       $corder = -1;
   }

#   if( ($strand == 1 && $contigb->orientation == 1) || ($strand == -1 && $contigb->orientation == -1) ) {
   if( $corder == 1 ) {
       $start = $contigb->offset - ($Bio::EnsEMBL::DB::CloneI::CONTIG_SPACING/2);
       $end = $start + $number -1;
       $isc = 0;
   } else {
       $start = $contigb->offset + $contigb->length + ($Bio::EnsEMBL::DB::CloneI::CONTIG_SPACING/2);
       $end = $start + $number -1;
       $isc = 1;
   }
   
   my $ret;
   if( $isc == 0 ) {
       $ret = "$start..$end";
   } else {
       $ret = "complement($start..$end)";
   }

   print STDERR "Generating $ret\n";

   return $ret;
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



















