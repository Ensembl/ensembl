
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

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  vsns-bcd-perl@lists.uni-bielefeld.de          - General discussion
  vsns-bcd-perl-guts@lists.uni-bielefeld.de     - Technically-oriented discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

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
  
  # build the map of how contigs go onto the vc coorindates
  $self->_build_contig_map($focus,$focusposition,$ori,$leftsize,$rightsize);

# set stuff in self from @args
  return $make; # success - we hope!
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

   $self->{'rightmostcontig'} = $current_contig->id();
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

