
#
# BioPerl module for Bio::EnsEMBL::DB::VirtualContig
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DB::VirtualContig - A virtual contig implementation 

=head1 SYNOPSIS

  #get a virtualcontig somehow

  # usual contig methods applicable:

  @features = $virtualcontig->get_all_SimilarityFeatures();
  @genes    = $virtualcontig->get_all_Genes();
  $seq      = $virtualcontig->seq();

  # extend methods

  # makes a new virtualcontig 5000 base pairs to the 5'
  $newvirtualcontig = $virtualcontig->extend(-5000,-5000);

  # makes a virtualcontig of maximal size
  $newvirtualcontig = $virtualcontig->extend_maximally();
  
=head1 DESCRIPTION

A virtual contig gives a contig interface that is built up
of RawContigs, and in which the features/genes which come
off them are in a single coordinate system. (genes may have
exons that occur outside the virtual contig).

This implementation is of the VirtualContig interface but is
a puer-perl implementation that can sit ontop of any 
RawContigI compliant object. For that reason I have put it
in Bio::EnsEMBL::DB scope, indicating that other database
implementations can use this object if they so wish.


=head1 AUTHOR - Ewan Birney

Email birney@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DB::VirtualContig;
use vars qw($AUTOLOAD @ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object

use Bio::Root::Object;
use Bio::EnsEMBL::DB::VirtualContigI;

my $VC_UNIQUE_NUMBER = 0;


@ISA = qw(Bio::Root::Object Bio::EnsEMBL::DB::VirtualContigI);

# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize(@args);

  my ($focus,$focusposition,$ori,$leftsize,$rightsize) = $self->_rearrange([qw( FOCUS FOCUSPOSITION ORI LEFT RIGHT)],@args);

  if( !defined $focus || !defined $focusposition || !defined $ori || !defined $leftsize || !defined $rightsize ) {
      $self->throw("Have to provide all arguments to virtualcontig, focus, focusposition, ori, left, right");
  }

  # set up hashes for the map
  $self->{'start'} = {};
  $self->{'startincontig'} = {};
  $self->{'contigori'} = {};
  
  # this actually stores the contig we are using
  $self->{'contighash'} = {};

  # this is for cache's of sequence features if/when we want them
  $self->{'_sf_cache'} = {};
  
  # build the map of how contigs go onto the vc coorindates
  $self->_build_contig_map($focus,$focusposition,$ori,$leftsize,$rightsize);

  $self->_unique_number($VC_UNIQUE_NUMBER++);

# set stuff in self from @args
  return $make; # success - we hope!
}


=head1 Implementations for the ContigI functions

These functions are to implement the ContigI interface


=head2 seq

 Title   : seq
 Usage   : $seq = $contig->seq();
 Function: Gets a Bio::PrimarySeqI object out from the contig
 Example :
 Returns : Bio::PrimarySeqI object
 Args    :


=cut

sub seq{
   my ($self) = @_;

   my $seq = $self->_seq_cache();
   if( defined $seq ) {
       return $seq;
   }

   # we have to move across the map, picking up the sequences,
   # truncating them and then adding them into the final product.

   my @contig_id = sort { $self->{'start'}->{$a} <=> $self->{'start'}->{$b} } keys %{$self->{'start'}};
   my $seq_string;
   foreach my $cid ( @contig_id ) {
       my $c = $self->{'contighash'}->{$cid};
       my $tseq = $c->seq();

       my $trunc;

       if( $self->{'contigori'}->{$c} == 1 ) {
	   $trunc = $tseq->subseq($self->{'start'}->{$cid},$c->golden_right);
       } else {
	   $trunc = $tseq->trunc($c->golden_left,$self->{'start'}->{$cid})->revcom->seq;
       }

       $seq_string .= $trunc;
   }

   $seq = Bio::PrimarySeq->new( -id => "virtual_contig_".$self->_unqiue_number,
				-seq => $seq_string,
				-moltype => 'dna'
				);


   $self->_seq_cache($seq);
  
   return $seq;
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

    return "virtual_contig_".$self->_unique_number;
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
   my @out;
   push(@out,$self->get_all_SimilarityFeatures());
   push(@out,$self->get_all_RepeatFeatures());
  # push(@out,$self->

   return @out;

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
   
   return $self->_get_all_SeqFeatures_type('similarity');

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
   
   return $self->_get_all_SeqFeatures_type('repeat');


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

   $self->throw("get_all_Genes on virtual contigs not implemented yet!");

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

   return $self->_left_size + $self->_right_size +1;

}




=head2 _build_contig_map

 Title   : _build_contig_map
 Usage   : Internal function for building the map
           of rawcontigs onto the virtual contig positions.

           To do this we need contig ids mapped to start positions
           in the vc, and the orientation of the contigs.
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _build_contig_map{
   my ($self,$focus,$focusposition,$ori,$left,$right) = @_;

   # we first need to walk down contigs going left
   # so we can figure out the start position (contig-wise)

   # initialisation - find the correct end of the focus contig
   my ($current_left_size,$current_orientation,$current_contig,$overlap);

   if( $ori == 1 ) {
       $current_left_size = $focusposition;
       $current_orientation = 1;
   } else {
       $current_left_size = $focus->length - $focusposition;
       $current_orientation = -1;
   }
   $current_contig = $focus;

   while( $current_left_size < $left ) {
       print STDERR "Looking at ",$current_contig->id," with $current_left_size\n";

       if( $current_orientation == 1 ) {
	   
	   # go left wrt to the contig.
	   $overlap = $current_contig->get_left_overlap();
	   # if there is no left overlap, trim left to this size
	   # as this means we have run out of contigs

	   print STDERR "Gone left\n";

	   if( !defined $overlap ) {
	       $left = $current_left_size;
	       print STDERR "getting out - no overlap\n";
	       last;
	   }
	   # add to total, move on the contigs
	   $self->warn("Not coping with non-overlapping, sized gaps");
	   $current_left_size += $overlap->sister->golden_length;
	   $current_contig = $overlap->sister();

	   if( $overlap->sister_polarity == 1) {
	       $current_orientation = 1;
	   } else {
	       $current_orientation = -1;
	   }
       } else {
	   # go right wrt to the contig.
	   $overlap = $current_contig->get_right_overlap();
	   # if there is no left overlap, trim left to this size
	   # as this means we have run out of contigs
	   if( !defined $overlap ) {
	       $left = $current_left_size;
	       last;
	   }

	   # add to total, move on the contigs
	   $self->warn("Not coping with non-overlapping, sized gaps");
	   $current_left_size += $overlap->sister->golden_length;
	   $current_contig = $overlap->sister();

	   if( $overlap->sister_polarity == 1) {
	       $current_orientation = -1;
	   } else {
	       $current_orientation = 1;
	   }
       }
   }



   # now $current_contig is the left most contig in this set, with
   # its orientation set and ready to rock... ;)

   my $total = $left + $right;

   print STDERR "leftmost contig is ",$current_contig->id,"with $total to account for\n";
   $self->{'leftmostcontig_id'} = $current_contig->id;

   # the first contig will need to be trimmed at a certain point
   my $startpos;
   if( $current_orientation == 1 ) {
       $startpos = $current_contig->golden_start + ($current_left_size - $left);
   } else {
       $startpos = $current_contig->golden_end   - ($current_left_size - $left);
   }

   print STDERR "Leftmost contig has $startpos and $current_orientation $left vs $current_left_size\n";

   $self->{'start'}->{$current_contig->id} = 1;
   $self->{'startincontig'}->{$current_contig->id} = $startpos;
   $self->{'contigori'}->{$current_contig->id} = $current_orientation;
   $self->{'contighash'}->{$current_contig->id} = $current_contig;

   
   my $current_length;

   if( $current_orientation == 1 ) {
       $current_length = $current_contig->golden_end - $startpos;
   } else {
       $current_length = $startpos - $current_contig->golden_start;
   }
   print STDERR "current length before we get into this is $current_length\n";


   while( $current_length < $total ) {
       print STDERR "In building actually got $current_length towards $total\n";
 
       # move onto the next contig.

       if( $current_orientation == 1 ) {
	   # go right wrt to the contig.
	   print STDERR "Going right\n";

	   $overlap = $current_contig->get_right_overlap();
	 
	   # if there is no right overlap, trim right to this size
	   # as this means we have run out of contigs
	   if( !defined $overlap ) {
	       print STDERR "Out of contigs!\n";

	       $right = $current_length - $left;
	       last;
	   }

	   # add to total, move on the contigs
	   $self->warn("Not coping with non-overlapping, sized gaps");

	   $current_contig = $overlap->sister();
	   $self->{'contighash'}->{$current_contig->id} = $current_contig;

	   if( $overlap->sister_polarity == 1) {
	       $current_orientation = 1;
	   } else {
	       $current_orientation = -1;
	   }

	   $self->{'start'}->{$current_contig->id} = $current_length; # maybe +1 (?)
	   if( $current_orientation == 1 ) {
	       $self->{'startincontig'}->{$current_contig->id} = $current_contig->golden_start;
	   } else {
	       $self->{'startincontig'}->{$current_contig->id} = $current_contig->golden_end;
	   }

	   $self->{'contigori'}->{$current_contig->id} = $current_orientation;

	   # up the length
	   $current_length += $overlap->sister->golden_length;
       } else {
	   # go left wrt to the contig
	   print STDERR "Going left\n";

	   $overlap = $current_contig->get_left_overlap();
	 
	   # if there is no left overlap, trim right to this size
	   # as this means we have run out of contigs
	   if( !defined $overlap ) {
	       $right = $current_length - $left;
	       last;
	   }

	   # add to total, move on the contigs
	   $self->warn("Not coping with non-overlapping, sized gaps");

	   $current_contig = $overlap->sister();
	   $self->{'contighash'}->{$current_contig->id} = $current_contig;

	   if( $overlap->sister_polarity == 1) {
	       $current_orientation = -1;
	   } else {
	       $current_orientation = 1;
	   }

	   $self->{'start'}->{$current_contig->id} = $current_length; # maybe +1 (?)
	   if( $current_orientation == 1 ) {
	       $self->{'startincontig'}->{$current_contig->id} = $current_contig->golden_start;
	   } else {
	       $self->{'startincontig'}->{$current_contig->id} = $current_contig->golden_end;
	   }

	   $self->{'contigori'}->{$current_contig->id} = $current_orientation;

	   # up the length
	   $current_length += $overlap->sister->golden_length;
       }
   }

   # $right might have been modified during the walk

   $total = $left + $right;

   # need to store end point for last contig

   $self->{'rightmostcontig_id'} = $current_contig->id();
   if( $current_orientation == 1 ) {
       $self->{'rightmostend'}    = $current_contig->golden_end - ($total - $current_length);
   } else {
       $self->{'rightmostend'}    = $current_contig->golden_start + ($total - $current_length);
   }
   
   # put away the focus/size info etc

   $self->_focus($focus);
   $self->_focus_position($focusposition);
   $self->_focus_orientation($ori);
   $self->_left_size($left);
   $self->_right_size($right);


   # ready to rock and roll. Woo-Hoo!

}





=head2 _get_all_SeqFeatures_type

 Title   : _get_all_SeqFeatures_type
 Usage   : Internal function which encapsulates getting
           features of a particular type and returning
           them in the VC coordinates, optionally cache'ing
           them.
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _get_all_SeqFeatures_type{
   my ($self,$type) = @_;

   if( $self->_cache_seqfeatures() && $self->_has_cached_type($type) ) {
       return $self->_get_cache($type);
   }

   # ok - build the sequence feature list...

   my $sf;
   if( $self->_cache_seqfeatures() ) {
       $sf = $self->_make_cache($type);
   } else {
       $sf = []; # will be destroyed when drops out of scope
   }


   foreach my $c ( values %{$self->{'_contighash'}} ) {
       if( $type eq 'repeat' ) {
	   push(@$sf,$c->get_all_RepeatFeatures());
       } elsif ( $type eq 'similarity' ) {
	   push(@$sf,$c->get_all_SimilarityFeatures());
       } elsif ( $type eq 'prediction' ) {
	   push(@$sf,$c->get_all_PredictionFeatures());
       } else {
	   $self->throw("Type $type not recognised");
       }
   }


   foreach $sf ( @$sf ) {
       $self->_convert_seqfeature_to_vc_coords($sf);
   }

   return @$sf;
}


=head2 _convert_seqfeature_to_vc_coords

 Title   : _convert_seqfeature_to_vc_coords
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _convert_seqfeature_to_vc_coords{
   my ($self,$sf) = @_;

   my $cid = $sf->seqname();
   if( !exists $self->{'contighash'}->{$cid} ) {
       $self->throw("Attempting to map a sequence feature with $cid on a virtual contig with no $cid");
   }

   if( $self->{'contigori'}->{$cid} == 1 ) {
       # ok - forward with respect to vc. Only need to add offset
       my $offset = $self->{'start'}->{$cid} - $self->{'startincontig'}->{$cid};
       $sf->start($sf->start + $offset);
       $sf->end($sf->start + $offset);
       # strand stays the same
   } else {
       my $offset = $self->{'start'}->{$cid} + $self->{'startincontig'}->{$cid};
       # flip strand
       $sf->strand($sf->strand * -1);
       
       # yup. A number of different off-by-one errors possible here

       my $tstart = $sf->start;
       my $tend   = $sf->end;

       $sf->start($offset - $tend);
       $sf->end($offset - $tstart);

   }

}


=head2 _dump_map

 Title   : _dump_map
 Usage   : Produces a dumped map for debugging purposes
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _dump_map{
   my ($self,$fh) = @_;

   ! defined $fh && do { $fh = \*STDERR};

   my @ids = keys %{$self->{'contighash'}};
   @ids = sort { $self->{'start'}->{$a} <=> $self->{'start'}->{$b} } @ids;

   foreach my $id ( @ids ) {
       print $fh "Contig $id starts:",$self->{'start'}->{$id}," start in contig ",$self->{'startincontig'}->{$id}," orientation ",$self->{'contigori'}->{$id},"\n";
   }
}


=head2 _focus

 Title   : _focus
 Usage   : $obj->_focus($newval)
 Function: 
 Example : 
 Returns : value of _focus
 Args    : newvalue (optional)


=cut

sub _focus{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'_focus'} = $value;
    }
    return $obj->{'_focus'};

}

=head2 _focus_position

 Title   : _focus_position
 Usage   : $obj->_focus_position($newval)
 Function: 
 Example : 
 Returns : value of _focus_position
 Args    : newvalue (optional)


=cut

sub _focus_position{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'_focus_position'} = $value;
    }
    return $obj->{'_focus_position'};

}

=head2 _focus_orientation

 Title   : _focus_orientation
 Usage   : $obj->_focus_orientation($newval)
 Function: 
 Example : 
 Returns : value of _focus_orientation
 Args    : newvalue (optional)


=cut

sub _focus_orientation{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'_focus_orientation'} = $value;
    }
    return $obj->{'_focus_orientation'};

}

=head2 _left_size

 Title   : _left_size
 Usage   : $obj->_left_size($newval)
 Function: 
 Example : 
 Returns : value of _left_size
 Args    : newvalue (optional)


=cut

sub _left_size{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'_left_size'} = $value;
    }
    return $obj->{'_left_size'};

}
=head2 _right_size

 Title   : _right_size
 Usage   : $obj->_right_size($newval)
 Function: 
 Example : 
 Returns : value of _right_size
 Args    : newvalue (optional)


=cut

sub _right_size{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'_right_size'} = $value;
    }
    return $obj->{'_right_size'};

}

=head2 _cache_seqfeatures

 Title   : _cache_seqfeatures
 Usage   : $obj->_cache_seqfeatures($newval)
 Function: 
 Returns : value of _cache_seqfeatures
 Args    : newvalue (optional)


=cut

sub _cache_seqfeatures{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_cache_seqfeatures'} = $value;
    }
    return $obj->{'_cache_seqfeatures'};

}

=head2 _has_cached_type

 Title   : _has_cached_type
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _has_cached_type{
   my ($self,$type) = @_;

   if ( exists $self->{'_sf_cache'}->{$type} ) {
       return 1;
   } else {
       return 0;
   }
}

=head2 _make_cache

 Title   : _make_cache
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _make_cache{
   my ($self,$type) = @_;

   if( $self->_has_cached_type($type) == 1) {
       $self->throw("Already got a cache for $type! Error in logic here");
   }

   $self->{'_sf_cache'}->{$type} = [];

   return $self->{'_sf_cache'}->{$type};
}

=head2 _get_cache

 Title   : _get_cache
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _get_cache{
   my ($self,$type) = @_;

   return $self->{'_sf_cache'}->{$type};
   
}

=head2 _unique_number

 Title   : _unique_number
 Usage   : $obj->_unique_number($newval)
 Function: 
 Returns : value of _unique_number
 Args    : newvalue (optional)


=cut

sub _unique_number{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_unique_number'} = $value;
    }
    return $obj->{'_unique_number'};

}

=head2 _seq_cache

 Title   : _seq_cache
 Usage   : $obj->_seq_cache($newval)
 Function: 
 Returns : value of _seq_cache
 Args    : newvalue (optional)


=cut

sub _seq_cache{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_seq_cache'} = $value;
    }
    return $obj->{'_seq_cache'};

}



1;
