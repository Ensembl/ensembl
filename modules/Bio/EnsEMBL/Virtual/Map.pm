
#
# Ensembl module for Bio::EnsEMBL::Virtual::Map
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Virtual::Map - Map of MapContigs which define a VirtualContig

=head1 SYNOPSIS

    # Virtual::Map objects are internal to VirtualContigs


=head1 DESCRIPTION

This object is basically a hash of MapContigs which make a virtual
contig  (MapContigs are just helper objects, see there).

=head1 AUTHOR - Ewan Birney

This modules is part of the Ensembl project http://www.ensembl.org

Email birney@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Virtual::Map;
use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Virtual::MapContig;

# Object preamble - inherits from Bio::Root::RootI

use Bio::Root::RootI;

@ISA = qw(Bio::Root::RootI);

# new() is written here 

sub new {
  my($class,@args) = @_;

  my $self = {};
  bless $self,$class;
  $self->{'_contig_map'} = {};
      
# set stuff in self from @args
  return $self;
}




=head2 build_map

 Title   : build_map
 Usage   : $map->build_map($rawcontig,$focusposition,$ori,$left,$right)
 Function: constructs a Map of the RawContigs in the map. Is (only?) called
   as one of the last step of the Virtual::Contig::new()
 Example :
 Returns : 
 Args    :


=cut

sub build_map {
   my ($self,$rawcontig,$focusposition,$ori,$left,$right) = @_;


   # pseudo-code
   # 1) calculate first offset
   # 2) walk left until left end condition met
   # 3) find last contig. Add to map
   # 4) walk right, added contigs to the map


   my $current_contig = $rawcontig;
   my $current_ori    = $ori;
   my $current_left_size = 0;
   

   if( $focusposition < $rawcontig->golden_start || $focusposition > $rawcontig->golden_end ) {
       $self->throw("Focus position is outside golden/start end on rawcontig. Ugh!");
   }


   if( $ori == 1 ) {
       $current_left_size = $focusposition - $rawcontig->golden_start;
   } else {
       $current_left_size = $rawcontig->golden_end - $focusposition;
   }

   my $left_overhang_size =0;

   # left walk
   while( 1 ) {
       print STDERR "Walking with $current_left_size to go compared to $left\n";

       
       if( $current_left_size >  $left ) {
	   $self->throw("Bad internal error. Did not terminate left walk on an explicit conition");
       }

       # go left
       my ($nextcontig,$start_in_contig,$contig_ori,$gap_distance) = 
	   $self->_go_left($current_contig,$current_ori);
       
       
       
       # if the gap pushes us over the total, we have a left overhang
       $current_left_size += $gap_distance;
       if( $current_left_size >= $left ) {
	   $self->left_overhang(1);
	   $left_overhang_size = $current_left_size + $gap_distance - $left;
	   ####
	   last;
       }
       # otherwise we want to include this contig
       
       # add gap distance to current_size. Remember, this could be 0
       $current_left_size += $gap_distance;
       
       if( $current_left_size + $nextcontig->golden_length < $left ) {
	   # add golden length distance
	   $current_left_size += $nextcontig->golden_length;
	   $current_contig = $nextcontig;
	   $current_ori    = $contig_ori;
	   next; # back to while(1)
       } else {
	   # this is the leftmost contig
	   $current_contig = $nextcontig;
	   
	   $self->left_overhang(0);
	   
       }

   }


   my $total = $left+$right;
   my $current_size = $left_overhang_size;

   # put in this contig as the first contig. 

   # main build

   while( 1 ) {

       # paranoia check - we should always exit on a defined condition
       if( $current_size >= $total ) {
	   $self->throw("Bad internal error: in right walk on vc build, got over distance without triggering a defined end condition");
       }

       # go right
       my ($nextcontig,$start_in_contig,$contig_ori,$gap_distance) = 
	   $self->_go_right($current_contig,$current_ori);
       
       # if the gap pushes us over the total, we have a right overhang
       $current_size += $gap_distance;
       if( $current_size >= $total ) {
	   $self->right_overhang(1);
	   last;
       }
       # otherwise we want to include this contig
       
       # add gap distance to current_size. Remember, this could be 0
       $current_size += $gap_distance;

       if( $current_size + $nextcontig->golden_length < $total ) {
	   # we want to include the entire contig
	   $self->create_MapContig($nextcontig,
				   $current_size+1,
				   $current_size+$nextcontig->golden_length,
				   $start_in_contig,
				   $contig_ori);
	   # add golden length distance
	   $current_size += $nextcontig->golden_length;
	   $current_contig = $nextcontig;
	   $current_ori    = $contig_ori;
	   next; # back to while(1)
       } else {
	   # this is the rightmost contig
	   my $length = $total - $current_size;

	   my $start;
	   if( $contig_ori == 1 ) {
	       $start = $start_in_contig;
	   } else {
	       # start has to move to right end - length
	       $start = $start_in_contig + $nextcontig->golden_length - $length;
	   }

	   $self->create_MapContig($nextcontig,
				   $current_size+1,
				   $current_size+$length,
				   $start,
				   $contig_ori);

	   $self->right_overhang(0);
	   last;
       }
   }
	   
}


=head2 create_MapContig

 Title   : create_MapContig
 Usage   :
 Function: creates a MapContig and adds it into own VirtualMap.
 Example :
 Returns : 
 Args    : RawContig, ChrStart, ChrEnd, RawStart, Orientation


=cut

sub create_MapContig {
   my ($self,$rawcontig,$start,$end,$start_in_rawcontig,$orientation) = @_;

   if( !defined $orientation || $start > $end ||
       !ref $rawcontig || !$rawcontig->isa('Bio::EnsEMBL::DB::ContigI') ) {
       $self->throw("Invalid arguments passed into create_MapContig ($rawcontig)");
   }

   my $out = Bio::EnsEMBL::Virtual::MapContig->new(
				       -rawcontig => $rawcontig,
				       -start => $start,
				       -end => $end,
				       -rawcontig_start => $start_in_rawcontig,
				       -orientation => $orientation
				       );

   $self->_add_MapContig($out);

   return $out;

}

=head2 get_MapContig_by_id

 Title   : get_MapContig_by_id
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_MapContig_by_id {
   my ($self,$name) = @_;

   return $self->{'_contig_map'}->{$name};
}

sub get_MapContig {
   my ($self,$name) = @_;

   my ($p,$f,$l) = caller;
   $self->warn("$f:$l:get_MapContig a deprecated method. use get_MapContig_by_id instead");
   return $self->get_MapContig_by_id($name);
}


=head2 each_MapContig

 Title   : each_MapContig
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub each_MapContig{
   my ($self,$reverse) = @_;

   my @mapcontigs = values %{$self->{'_contig_map'}};

   if ($reverse) {
       @mapcontigs = sort { $b->start <=> $a->start} @mapcontigs;
   }
   else {
       @mapcontigs = sort { $a->start <=> $b->start} @mapcontigs;
   }

   return (@mapcontigs);
}

sub get_all_MapContigs {
    my ($self) = @_;

    my ($p,$f,$l) = caller;
    $self->warn("$f:$l:get_all_MapContigs a deprecated method. use each_MapContig instead");

    return $self->each_MapContig;
}


sub get_all_RawContigs {
    my ($self) = @_;

    my @contigs;

    foreach my $contig ($self->each_MapContig) {
	push(@contigs,$contig->contig);
    }

    return @contigs;
}

=head2 raw_contig_position

 Title   : raw_contig_position
 Usage   : my ($map_contig,$rc_position,$rc_strand) = $vmap->raw_contig_position($vc_pos,$vc_strand)
 Function: Maps a VirtualContig position to the RawContig Position
 Returns : Bio::EnsEMBL::Virtual::MapContig object, 
           position (int), strand (int)
 Args    : position (int), strand (int)


=cut

sub raw_contig_position {
    my ($self, $vcpos, $vcstrand)=@_;
 
    my $rc;
    my $rc_pos;
    my $rc_strand;
    
    #my $length=$self->length;
    #if ($vcpos >$length) {
    #	$self->throw("Asked to map vc position outside vc coordinates!\n");
    #}
    #print STDERR "Looking for $vcpos....\n";
    #Go through all Contigs and find out where vcpos lies

    foreach my $mc ($self->each_MapContig) {
	
	#If we are still on a RawContig which finished vcpos, 
        #move to next contig
	if ($mc->end < $vcpos) {
	    next;
	}
	
	#If vcpos is within the start and enf of this Contig, we found it!
	#And we get out of the loop...
	if (($vcpos >= $mc->start)&&($vcpos <= $mc->end)) {
	    #print STDERR "Found contig!\n ".$mc->contig->id." with start ".$mc->start." and end ".$mc->end."\n"; 
	    $rc=$mc->contig;
	    if ($mc->orientation == 1) {
		# the contig starts at startin
		$rc_pos = $mc->rawcontig_start + ($vcpos - $mc->start);
	    }
	    else {
		$rc_pos=$mc->rawcontig_end - ( $vcpos-$mc->start );
	    }
	    
	    #If strand passed to the method, sort out the strand
	    if ($vcstrand) {
		if ($vcstrand == 1) {
		    $rc_strand = $mc->orientation;
		}
		else {
		    $rc_strand = -$mc->orientation;
		}
	    }
	    last;
	}

	#If we are not out of the loop at this stage it means that
	#our Contig lies in a gap
	if ($mc->start > $vcpos) {
	    $rc='N';
	    $rc_pos=-1;
	    $rc_strand=-2;
	    last;
	}
    }
    
    
    $vcstrand && return $rc,$rc_pos,$rc_strand;
    return $rc,$rc_pos;
}

sub vcpos_to_rcpos {
    my ($self, $vcpos, $vcstrand)=@_;
    my ($p,$f,$l) = caller;
    $self->warn("$f:$l:vcpos_to_rcpos a deprecated method. use raw_contig_position instead");
    return $self->raw_contig_position($vcpos,$vcstrand);
}

    

=head2 length

 Title   : length
 Usage   : $obj->length($newval)
 Function: 
 Example : 
 Returns : value of length
 Args    : newvalue (optional)


=cut

sub length{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'length'} = $value;
    }
    return $obj->{'length'};

}



=head2 right_overhang

 Title   : right_overhang
 Usage   : $obj->right_overhang($newval)
 Function: 
 Example : 
 Returns : value of right_overhang
 Args    : newvalue (optional)


=cut

sub right_overhang{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'right_overhang'} = $value;
    }
    return $obj->{'right_overhang'};

}

=head2 left_overhang

 Title   : left_overhang
 Usage   : $obj->left_overhang($newval)
 Function: 
 Example : 
 Returns : value of left_overhang
 Args    : newvalue (optional)


=cut

sub left_overhang{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'left_overhang'} = $value;
    }
    return $obj->{'left_overhang'};

}


=head2 _add_MapContig

 Title   : _add_MapContig
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _add_MapContig{
   my ($self,$mc) = @_;

   $self->{'_contig_map'}->{$mc->contig->id} = $mc;
}

=head2 _go_left

 Title   : _go_left
 Usage   : ($nextcontig,$nextcontigpos,$orientation,$distance) = $mc->_go_left($current_contig,$ori)
 Function: A helper function for traversing rawcontig chains. Should only
           be used by people who know what they are doing
 Example :
 Returns : 
 Args    :


=cut

sub _go_left{
   my ($self,$current,$ori) = @_;

   if( $ori == 1 ) {
       my $co = $current->get_left_overlap();
       my $start_in_contig;
       return ($co->sister,$co->sister->golden_start,$co->sister_polarity,$co->distance);
   } else {
       my $co = $current->get_right_overlap();
       my $pol = $co->sister_polarity * -1;
       return ($co->sister,$co->sister->golden_start,$pol,$co->distance);
   }

}

=head2 _go_right

 Title   : _go_right
 Usage   : ($nextcontig,$nextcontigpos,$orientation,$distance) = $mc->_go_right($current_contig,$ori)
 Function: A helper function for traversing rawcontig chains. Should only
           be used by people who know what they are doing
 Example :
 Returns : 
 Args    :


=cut

sub _go_right{

   my ($self,$current,$ori) = @_;

   if( $ori == 1 ) {
       my $co = $current->get_right_overlap();
       return ($co->sister,$co->sister->golden_start,$co->sister_polarity,$co->distance);
   } else {
       my $co = $current->get_left_overlap();
       my $pol = $co->sister_polarity * -1;
       return ($co->sister,$co->sister->golden_start,$pol,$co->distance);
   }

}


1;



