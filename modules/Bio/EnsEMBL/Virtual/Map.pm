
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

   if (!defined $self->{'_contig_map'}->{$name}) {
       $self->throw("Could not find Map Contig for contig $name\n");
   }
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
 Function: Maps a VirtualContig position to a RawContig + RawContig position

 Returns : The underlying RawContig and a position on it (in RC coords),
           and optionally the RC strandedness
 Args   : position on VirtualContig (in VC coords), and optionally
          VirtualContig strand.
=cut

sub raw_contig_position {
# PL: belongs in Contig? 
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
	    $rc='gapcontig';
	    $rc_pos= $vcpos;
	    $rc_strand= $vcstrand;
	    last;
	}
    }
    
    $vcstrand && return $rc,$rc_pos,$rc_strand;
    return $rc,$rc_pos;
}

=head2 raw_contig_interval

 Title   : raw_contig_interval
 Usage   : my ($map_contig,$rc_position,$rc_strand) = $vmap->raw_contig_position($vc_pos,$vc_strand)
 Function: Maps a VirtualContig position to a RawContig + RawContig position

 Returns : The underlying RawContig and a position on it (in RC coords),
           and optionally the RC strandedness
 Args   : position on VirtualContig (in VC coords), and optionally
          VirtualContig strand.
=cut

sub raw_contig_interval {

    my ($self, $vc_start, $vc_end, $vc_strand)=@_;
 
    if( ! defined $vc_strand ) {
      $self->throw( "Wrong number of arguments. " );
    }

    my $rc;
    my $rc_pos;
    my $rc_strand;
    

    my @result;
    my $last_mc;

    foreach my $mc ($self->each_MapContig) {
	
      my ( $raw_id, $raw_start, $raw_end, $raw_strand  );
	#If we are still on a RawContig which finished vcpos, 
        #move to next contig
	if ($mc->end < $vc_start) {
	    next;
	}
	
      if ( $mc->start > $vc_end ) {
	last;
      }

      if ( $vc_start < $mc->start ) {
	# ups, gap detected
	push( @result, { "gap_start" => $vc_start,
			 "gap_end" => $mc->start-1 } );
	$vc_start = $mc->start;
      }
      # VC start is somewhere in mapcontig region
      if ($mc->orientation == 1) {
	# 
	$raw_start = $mc->rawcontig_start + ($vc_start - $mc->start);
      } else {
	$raw_end = $mc->rawcontig_end - ( $vc_start - $mc->start );
      }

      if( $vc_end > $mc->end ) {
	if( $mc->orientation == 1 ) {
	  $raw_end = $mc->rawcontig_end;
	} else {
	  $raw_start = $mc->rawcontig_start;
	}
      } else {
	if ($mc->orientation == 1) {
	  $raw_end = $mc->rawcontig_start + ($vc_end - $mc->start);
	} else {
	  $raw_start = $mc->rawcontig_end - ( $vc_end - $mc->start );
	}
      }
      push( @result, { "raw_contig_id" => $mc->contig->internal_id, 
		      "raw_start" => $raw_start, 
		       "raw_end" => $raw_end, 
		       "raw_strand" => $mc->orientation * $vc_strand } );
	$last_mc = $mc;
       	$vc_start = $mc->end+1;
    }

    if( !defined($last_mc) ) {
      push( @result, { "gap_start" => $vc_start,
		       "gap_end" => $vc_end } ); 
    } elsif( $last_mc->end < $vc_end ) {
      # there was a gap at the end
      push( @result, { "gap_start" => $last_mc->end+1,
		       "gap_end" => $vc_end } );
    }
    
    if( $vc_strand == -1 ) {
      @result = reverse( @result );
    }

    return @result;
}


sub vcpos_to_rcpos {
#PL: belongs in Contig.pm ? 
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



