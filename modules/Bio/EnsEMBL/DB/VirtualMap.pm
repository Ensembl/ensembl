#
# BioPerl module for VirtualMap
#
# Cared for by Elia Stupka <elia@ebi.ac.uk>
#
# Copyright Elia Stupka
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

VirtualMap - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DB::VirtualMap;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;
use Bio::EnsEMBL::DB::MapContig;

@ISA = qw(Bio::Root::Object);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self) = @_;
  
  my $make = $self->SUPER::_initialize;
  $self->{'mapcontighash'}= {};

# set stuff in self from @args
 return $make; # success - we hope!
}

=head2 get_MapContig

 Title   : get_MapContig
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_MapContig{
   my ($self,$id) = @_;

   $id || $self->throw("Need to give an id to get a mapcontig!\n");
   
   if (my $mapcontig = $self->{'mapcontighash'}->{$id}) {
       return $mapcontig;
   }
   else {
       $self->throw("Could not find mapcontig $id");
   }
}

=head2 add_MapContig

 Title   : add_MapContig
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub create_MapContig{
   my ($self,$start,$start_in,$ori,$contig) = @_;

   $start || $self->throw("Need to give a start for the MapContig!");
   $start_in || $self->throw("Need to give a start_in for the MapContig!");
   $ori || $self->throw("Need to give an orientation for the MapContig!");

   if( ! $contig->isa("Bio::EnsEMBL::DB::RawContigI") ) {
       $self->throw("$contig is not a Bio::EnsEMBL::DB::RawContig!");
   }
   
   my $id=$contig->id;
   my $mapcontig=Bio::EnsEMBL::DB::MapContig->new( -contig =>$contig,
						   -ori => $ori,
						   -start => $start,
						   -startin => $start_in
						    );
   $self->{'mapcontighash'}->{$id}=$mapcontig;
}

=head2 get_all_MapContigs

 Title   : get_all_MapContigs
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_MapContigs{
   my ($self,$reverse) = @_;
   
   my @mapcontigs = values %{$self->{'mapcontighash'}};
   if ($reverse) {
       @mapcontigs = sort { $b->start <=> $a->start} @mapcontigs;
   }
   else {
       @mapcontigs = sort { $a->start <=> $b->start} @mapcontigs;
   }
   return (@mapcontigs);
}

=head2 get_all_RawContigs

 Title   : get_all_RawContigs
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_RawContigs{
   my ($self) = @_;
   
   my @contigs;

   foreach my $mc ($self->get_all_MapContigs){
       push @contigs,$mc->contig;
   }
   return (@contigs);
}

=head2 get_all_RawContig_ids

 Title   : get_all_Raw_Contig_ids
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub RawContig_ids{
   my ($self) = @_;
   
   return keys %{$self->{'mapcontighash'}};
}



=head2 build_clone_map

 Title   : build_clone_map
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub build_clone_map {
    my ($self,$clone) = @_;
    
    my $total_len   = 0;
    my $length      = 0;
    my $seen        = 0;
    my $middle      = 0;
    
    foreach my $contig ( $clone->get_all_Contigs ) {
	$self->create_MapContig($contig->embl_offset,1,1,$contig);
	
	$total_len = $contig->embl_offset + $contig->length;
	
	if( $total_len > $length ) {
	    $length = $total_len;
	}
	
	if( $seen == 0 ) {
	    $self->dbobj($contig->dbobj);
	    $seen = 1;
	}
	
    }
   
    # Tony: This vc made from a clone. Since it must have a left/right arm
    # we set the 'focus' to the middle.
    # The magic -1 avoids counting the focus base twice
    $middle = int($length)/2;
    $self->left_size($middle-1);
    $self->right_size($length-$middle);
    
    # Remember this vc contructed from a clone (rather than extending a 'seed' contig)
    $self->clone_map(1);
}


=head2 build_contig_map

 Title   : build_contig_map
 Usage   : Internal function for building the map
           of rawcontigs onto the virtual contig positions.

           To do this we need contig ids mapped to start positions
           in the vc, and the orientation of the contigs.
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub build_contig_map {
    my ($self,$focuscontig,$focusposition,$ori,$left,$right) = @_;
    
    # we first need to walk down contigs going left
    # so we can figure out the start position (contig-wise)
    # initialisation - find the correct end of the focus contig
    
    my ($current_left_size,$current_orientation,$current_contig,$overlap);
    $current_contig = $focuscontig;
    if( $focusposition < $current_contig->golden_start || $focusposition > $current_contig->golden_end ) {
	$self->throw("focus position is before or after golden region. Focus $focusposition [".$current_contig->golden_start.":".$current_contig->golden_end);
}
 
    if( $ori == 1 ) {
	$current_left_size   = $focusposition - $current_contig->golden_start;
	$current_orientation = 1;
    } else {
	$current_left_size   = $focuscontig->golden_end - $focusposition;
	$current_orientation = -1;
    }
    my %seen_hash;

    #print STDERR "Starting with $current_left_size and $left to build\n";

    GOING_LEFT :
    
	while( $current_left_size < $left ) {
	    #print STDERR "Current contig ".$current_contig->id." current left $current_left_size vs $left\n";

	    if( $seen_hash{$current_contig->id} ) {
		$self->throw("Bad internal error. Managed to loop back to the same contig in a virtualcontig walk. Something is inconsistent in the database. Id:".$current_contig->id);
	    } else {
		$seen_hash{$current_contig->id} = 1;
	    }
	
	    if( $current_orientation == 1 ) {
	      
		# go left wrt to the contig.
		$overlap = $current_contig->get_left_overlap();
		
		# if there is no left overlap, trim left to this size
		# as this means we have run out of contigs
	      
		if( !defined $overlap ) {
		    $left = $current_left_size;
                    #print STDERR "run out of contigs goin left\n";
		    last;
		}
		
		if( $overlap->distance == 1 ) {
		    $current_left_size += $overlap->sister->golden_length -1;
		} else {
		    $current_left_size += $overlap->distance;
		    if( $current_left_size > $left ) {
			# set the left overhang!
			$self->left_overhang($overlap->distance - ($current_left_size - $left));
			last GOING_LEFT;
		    }
		    
		    $current_left_size += $overlap->sister->golden_length;
		}
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
                    #print STDERR "run out of contigs going right\n";
		    last;
		}
		
		if( $overlap->distance == 1 ) {
		    $current_left_size += $overlap->sister->golden_length-1;
		} else {
		    $current_left_size += $overlap->distance;
		    if( $current_left_size > $left ) {
			# set the left overhang!
			$self->left_overhang($overlap->distance - ($current_left_size - $left));
			last GOING_LEFT;
		    }
		    $current_left_size += $overlap->sister->golden_length;
		}
		
		$current_contig = $overlap->sister();
		
		if( $overlap->sister_polarity == 1) {
		    $current_orientation = -1;
		} else {
		    $current_orientation = 1;
		}
	    }
            #print STDERR "current_left_size = $current_left_size\n";
	}
  
    # now $current_contig is the left most contig in this set, with
    # its orientation set and ready to rock... ;)
  
    my $total = $left + $right;
  
    # the first contig will need to be trimmed at a certain point
    my $startpos;
    
    my $current_length;

    if( $self->left_overhang() == 0 ) {
	if( $current_orientation == 1 ) {
	    $startpos = $current_contig->golden_start + ($current_left_size - $left);
	} else {
	    $startpos = $current_contig->golden_end   - ($current_left_size - $left);
	}
	# mysterious +1 to keep the overlapping base convention working.
	$self->create_MapContig(1,$startpos,$current_orientation,$current_contig);
	my $mc=$self->get_MapContig($current_contig->id);
	$mc->leftmost(1);
	
	if( $current_orientation == 1 ) {
	    $current_length = $current_contig->golden_end - $startpos +1;
	} else {
	    $current_length = $startpos - $current_contig->golden_start+1;
	}
    } else {
	# has an overhang - first contig offset into the system

	my $mc_start=$self->left_overhang+1;
	
	my $mc_startin;
	if( $current_orientation == 1 ) {
	    $mc_startin=$current_contig->golden_start;
	} else {
	    $mc_startin=$current_contig->golden_end;
	}

	$self->create_MapContig($mc_start,$mc_startin,$current_orientation,$current_contig);
	my $mc=$self->get_MapContig($current_contig->id);
	$mc->leftmost(1);
	
	$current_length = $self->left_overhang() + $current_contig->golden_length ;
    }
    
    # flush $seen_hash
    %seen_hash = ();
    
    while( $current_length < $total ) {
	#print STDERR "Looking on right move $current_length vs $total\n";

	# move onto the next contig.
	
	if( $seen_hash{$current_contig->id} ) {
	    $self->throw("Bad internal error. Managed to loop back to the same contig in a virtualcontig walk. Something is inconsistent in the database. Id:".$current_contig->id);
	} else {
	    $seen_hash{$current_contig->id} = 1;
	}
	
	if( $current_orientation == 1 ) {
	    # go right wrt to the contig.
	    
	    $overlap = $current_contig->get_right_overlap();
	    
	    # if there is no right overlap, trim right to this size
	    # as this means we have run out of contigs
	    if( !defined $overlap ) {
		#print STDERR "Found right end\n";
		$self->found_right_end(1);
		$right = $current_length - $left;
		last;
	    }
	    
	    # see whether the distance gives us an end condition, and a right_overhang
	    
	    if( $current_length + $overlap->distance > $total ) {
		#print STDERR "Found right overhang\n";
		# right overhang
		$self->right_overhang($total - $current_length);
		last;
	    }
	    
	    # add to total, move on the contigs
	    
	    $current_contig = $overlap->sister();
	   
	    if( $overlap->sister_polarity == 1) {
		$current_orientation = 1;
	    } else {
		$current_orientation = -1;
	    }
	    
	    # The +1, ++ and -- 's here are to handle the fact we want to produce 
            #abuting coordinate systems from overlapping switch points.
	    
	    my $mc_start;
	    if( $overlap->distance == 1 ) {
		$mc_start=$current_length +1;
	    } else {
		$mc_start=$current_length + $overlap->distance;
		$current_length += $overlap->distance;
	    }
	    
	    my $mc_startin;
	    if( $current_orientation == 1 ) {
		$mc_startin=$current_contig->golden_start;
		($overlap->distance == 1 ) && $mc_startin++; 
	    } else {
		$mc_startin=$current_contig->golden_end;
		( $overlap->distance == 1 ) && $mc_startin--;
	    }
	    
	    $self->create_MapContig($mc_start,$mc_startin,$current_orientation,$current_contig);
	    
	    # up the length
	    $current_length += $overlap->sister->golden_length -1;
	    #print STDERR "Returning with length $current_length\n";
	} else {
	    # go left wrt to the contig
	 	    
	    $overlap = $current_contig->get_left_overlap();

	    # if there is no left overlap, trim right to this size
	    # as this means we have run out of contigs
	    #IS THIS FINE?

	    if( !defined $overlap ) {
		#print STDERR "Found right end...going left\n";
		$self->found_right_end(1);
		$right = $current_length - $left;
		last;
	    }
	    
	    # see whether the distance gives us an end condition, and a right_overhang
	    
	    if( $current_length + $overlap->distance > $total ) {
		#print STDERR "Found right overhang\n";
		# right overhang
		$self->right_overhang($total - $current_length);
		last;
	    }

	    # add to total, move on the contigs
	    $current_contig = $overlap->sister();
	    
	    if( $overlap->sister_polarity == 1) {
		$current_orientation = -1;
	    } else {
		$current_orientation = 1;
	    }
	    
	    # The +1's here are to handle the fact we want to produce abutting
	    # coordinate systems from overlapping switch points.
	    my $mc_start;
	    if( $overlap->distance == 1 ) {
		$mc_start=$current_length +1;
	    } else {
		$mc_start=$current_length + $overlap->distance;
		$current_length += $overlap->distance;
	    }

	    my $mc_startin;
	    if( $current_orientation == 1 ) {
		$mc_startin=$current_contig->golden_start+1;
	    } else {
		$mc_startin=$current_contig->golden_end-1;
	    }
	    
	    $self->create_MapContig($mc_start,$mc_startin,$current_orientation,$current_contig);
	    
	    # up the length
	    $current_length += $overlap->sister->golden_length -1;
	    #print STDERR "Returning with length $current_length\n";
	}
    }
    
    # $right might have been modified during the walk
    
    $total = $left + $right;

    # need to store end point for last contig
    my $mc=$self->get_MapContig($current_contig->id);
   
    my $end;
    
    if( $self->right_overhang == 0 ) {
	if( $current_orientation == 1 ) {
	    $end= $current_contig->golden_end - ($current_length - $total);
	} else {
	    $end= $current_contig->golden_start + ($current_length - $total);
	}
    } else {
	if( $current_orientation == 1 ) {
	    $end= $current_contig->golden_end;
	} else {
	    $end= $current_contig->golden_start;
	}
    }
    $mc->rightmost_end($end);
   
    # put away the focus/size info etc

    $self->focus_contig($focuscontig);
    $self->focus_position($focusposition);
    $self->focus_orientation($ori);
    $self->left_size($left);
    $self->right_size($right);
    
    # ready to rock and roll. Woo-Hoo!
}

=head2 focus_contig

 Title   : focus_contig
 Usage   : $obj->focus_contig($newval)
 Function: 
 Example : 
 Returns : value of focus_contig
 Args    : newvalue (optional)


=cut

sub focus_contig {
    my ($obj,$value) = @_;
    if( defined $value) {
	$obj->{'_focus_contig'} = $value;
    }
    return $obj->{'_focus_contig'};
}

=head2 focus_position

 Title   : focus_position
 Usage   : $obj->focus_position($newval)
 Function: 
 Example : 
 Returns : value of focus_position
 Args    : newvalue (optional)


=cut

sub focus_position {
    my ($obj,$value) = @_;
    if( defined $value) {
	$obj->{'_focus_position'} = $value;
    }
    return $obj->{'_focus_position'};
}

=head2 focus_orientation

 Title   : focus_orientation
 Usage   : $obj->focus_orientation($newval)
 Function: 
 Example : 
 Returns : value of focus_orientation
 Args    : newvalue (optional)


=cut

sub focus_orientation {
    my ($obj,$value) = @_;
    if( defined $value) {
	$obj->{'_focus_orientation'} = $value;
    }
    return $obj->{'_focus_orientation'};
}

=head2 left_size

 Title   : left_size
 Usage   : $obj->left_size($newval)
 Function: 
 Example : 
 Returns : value of left_size
 Args    : newvalue (optional)


=cut

sub left_size {
    my ($obj,$value) = @_;
    if( defined $value) {
	$obj->{'_left_size'} = $value;
    }
    return $obj->{'_left_size'};
}

=head2 right_size

 Title   : right_size
 Usage   : $obj->right_size($newval)
 Function: 
 Example : 
 Returns : value of right_size
 Args    : newvalue (optional)


=cut

sub right_size {
    my ($obj,$value) = @_;
    if( defined $value) {
	$obj->{'_right_size'} = $value;
    }
    return $obj->{'_right_size'};
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
	$obj->{'_left_overhang'} = $value;
    }
    return $obj->{'_left_overhang'};   
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
	$obj->{'_right_overhang'} = $value;
    }
    return $obj->{'_right_overhang'};   
}

=head2 clone_map

 Title   : clone_map
 Usage   : $obj->clone_map($newval)
 Function: 
 Example : 
 Returns : value of clone_map
 Args    : newvalue (optional)


=cut

sub clone_map {
    my ($obj,$value) = @_;
    if( defined $value) {
	$obj->{'_clone_map'} = $value;
    }
    return $obj->{'_clone_map'};
}


=head2 dbobj

 Title   : dbobj
 Usage   : $obj->dbobj($newval)
 Function: 
 Example : 
 Returns : value of dbobj
 Args    : newvalue (optional)


=cut

sub dbobj{
    my ($obj,$value) = @_;
    
    if( defined $value) {
	$obj->{'dbobj'} = $value;
    }
    return $obj->{'dbobj'};   
}

=head2 found_right_end

 Title   : found_right_end
 Usage   : 
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub found_right_end {
    my ($self,$arg) = @_;
    
    if (defined($arg) && ($arg == 1 || $arg == 0)) {
	$self->{'_found_right_end'} = $arg;
    } elsif (defined($arg)) {
	$self->throw("Arg to found_right_end should be 0,1");
    }
    
    return $self->{'_found_right_end'};
}

=head2 found_left_end

 Title   : found_left_end
 Usage   : 
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub found_left_end {
    my ($self, $arg) = @_;
    
    if (defined($arg) && ($arg == 1 || $arg == 0)) {
	
	$self->{'_found_left_end'} = $arg;
    } elsif (defined($arg)) {
	$self->throw("Arg to found_left_end should be 0,1");
    }
    
    return $self->{'_found_left_end'};
}

=head2 vcpos_to_rcpos

 Title   : vcpos_to_rcpos
 Usage   : Deprecated: use raw_contig_position instead


=cut

sub vcpos_to_rcpos {
    my $self = shift;
    $self->warn("vcpos_to_rcpos: Deprecated name, use raw_contig_position instead\n");
    $self->raw_contig_position(@_);
}

=head2 raw_contig_position

 Title   : raw_contig_position
 Usage   : my ($map_contig,$rc_position,$rc_strand) = $vmap->raw_contig_position($vc_pos,$vc_strand)
 Function: Maps a VirtualContig position to the RawContig Position
 Returns : Bio::EnsEMBL::DB::MapContig object, 
           position (int), strand (int)
 Args    : position (int), strand (int)


=cut

sub raw_contig_position {
    my ($self, $vcpos, $vcstrand)=@_;
 
    my $rc;
    my $rc_pos;
    my $rc_strand;
    
    my $length=$self->left_size + $self->right_size;
    if ($vcpos >$length) {
	$self->throw("Asked to map vc position outside vc coordinates!\n");
    }
    #print STDERR "Looking for $vcpos....\n";
    #Go through all Contigs and find out where vcpos lies
    foreach my $mc ($self->get_all_MapContigs) {
	
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
		$rc_pos = $mc->start_in + ($vcpos - $mc->start);
		#$rc_pos=$vcpos-$mc->start+$mc->start_in;
	    }
	    else {
		# the contig starts at start in but reversed, hence the subtraction
		$rc_pos=$mc->start_in - ( $vcpos-$mc->start );
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

