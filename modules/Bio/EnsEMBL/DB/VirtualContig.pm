
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
 
  $vc = Bio::EnsEMBL::DB::VirtualContig->new( -focuscontig => $rawcontig,
					      -focusposition => 2,
					      -ori => 1,
					      -left => 5000,
					      -right => 5000
					      );

  # or
  $vc = Bio::EnsEMBL::DB::VirtualContig->new( -clone => $clone);

  # usual contig methods applicable:

  @features   = $virtualcontig->get_all_SimilarityFeatures();
  @genes      = $virtualcontig->get_all_Genes();
  $seq        = $virtualcontig->primary_seq();

  # extend methods

  # makes a new virtualcontig 5000 base pairs to the 5'
  $newvirtualcontig = $virtualcontig->extend(-5000,-5000);

  # makes a virtualcontig of maximal size
  $newvirtualcontig = $virtualcontig->extend_maximally();
  
=head1 DESCRIPTION

A virtual contig gives a contig interface that is built up
of RawContigs where the features/genes that come
off them are in a single coordinate system. (genes may have
exons that occur outside the virtual contig).

This implementation is of the VirtualContig interface but is
a pure-perl implementation that can sit atop any 
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

  $self->name("Virtual Contig Module"); # set the exception context (does this work?)

  my ($focuscontig,$focusposition,$ori,$leftsize,$rightsize,$clone) = $self->_rearrange([qw( FOCUSCONTIG FOCUSPOSITION ORI LEFT RIGHT CLONE)],@args);

  # set up hashes for the map
  $self->{'start'}         = {};
  $self->{'startincontig'} = {};
  $self->{'contigori'}     = {};
  $self->{'feature_skip'}  = {};
  # this actually stores the contig we are using
  $self->{'contighash'}    = {};

  # this is for cache's of sequence features if/when we want them
  $self->{'_sf_cache'}     = {};

  $self->_left_overhang(0);
  $self->_right_overhang(0);

  if( defined $clone && defined $focuscontig ){
      $self->throw("Build a virtual contig either with a clone or a focuscontig, but not with both");
  }
  
  if( defined $clone ) {
      $self->_build_clone_map($clone);
      $VC_UNIQUE_NUMBER = $clone->id;
      
  } else {
      if( !defined $focuscontig   || 
	  !defined $focusposition || 
	  !defined $ori           || 
	  !defined $leftsize      || 
	  !defined $rightsize ) {
	  $self->throw("Have to provide all arguments to virtualcontig \n" .
		       "(focuscontig, focusposition, ori, left and right)");
      }
      
      # build the map of how contigs go onto the vc coorindates
      $self->_build_contig_map($focuscontig,$focusposition,$ori,$leftsize,$rightsize);
      $self->dbobj($focuscontig->dbobj);
      $VC_UNIQUE_NUMBER = $focuscontig->id.".$focusposition.$ori.$leftsize.$rightsize";
  }
  
  $self->_unique_number($VC_UNIQUE_NUMBER);
  

# set stuff in self from @args
  return $make; # success - we hope!
}


=head1 Implementations for the ContigI functions

These functions are to implement the ContigI interface


=head2 extend

 Title   : extend
 Usage   : $new_vc = $vc->extend(100,100);
 Function: Make a new vc by extending an existing one
 Example :
 Returns : Bio::EnsEMBL::DB::VirtualContig
 Args    :


=cut

sub extend {
   my ($self, $left, $right) = @_;

   if( !ref $self || ! $self->isa('Bio::EnsEMBL::DB::VirtualContigI') ) {
       $self->throw("Can only extend a VirtualContigI, Bailing out...");
   }

   if (! defined $left || ! defined $right ){
       $self->throw("Must supply a left and right value when extending a VirtualContig");
   }
   
   print STDERR "Extending raw contig ".$self->_focus_contig->id." (ori = $self->_focus_orientation)\n";

   my $nvc = Bio::EnsEMBL::DB::VirtualContig->new( -focuscontig => $self->_focus_contig,
					       -focusposition   => $self->_focus_position,
					       -ori             => $self->_focus_orientation,
					       -left            => $self->_left_size - $left,
					       -right           => $self->_right_size + $right,
					       );


   my $id = join('.', ($nvc->_focus_contig->id, $nvc->_focus_position, $nvc->_focus_orientation, $nvc->_left_size, $nvc->_right_size));
   $nvc->_unique_number($id);
   
   return $nvc;
}


=head2 extend_maximally

 Title   : extend_maximally
 Usage   : $new_vc = $vc->extend_maximally();
 Function: Extends an existing vc as far as possible in both directions
 Example :
 Returns : Bio::EnsEMBL::DB::VirtualContig
 Args    :


=cut

sub extend_maximally {
   my ($self) = @_;

   if( !ref $self || ! $self->isa('Bio::EnsEMBL::DB::VirtualContigI') ) {
       $self->throw("Can only extend a VirtualContigI, Bailing out...");
   }
   # based on an original idea by Ewan Birney. ;)
   my $nvc = $self->extend(10000000000,10000000000);
   return $nvc;
}


=head2 extend_maximally_left

 Title   : extend_maximally_left
 Usage   : $new_vc = $vc->extend_maximally_left();
 Function: Extends an existing vc as far as possible to the left
 Example :
 Returns : Bio::EnsEMBL::DB::VirtualContig
 Args    :


=cut

sub extend_maximally_left {
   my ($self) = @_;

   if( !ref $self || ! $self->isa('Bio::EnsEMBL::DB::VirtualContigI') ) {
       $self->throw("Can only extend a VirtualContigI, Bailing out...");
   }
   # based on an original idea by Ewan Birney. ;)
   my $nvc = $self->extend(10000000000,0);
   return $nvc;
}


=head2 extend_maximally_right

 Title   : extend_maximally_right
 Usage   : $new_vc = $vc->extend_maximally_right();
 Function: Extends an existing vc as far as possible to the right
 Example :
 Returns : Bio::EnsEMBL::DB::VirtualContig
 Args    :


=cut

sub extend_maximally_right {
   my ($self) = @_;

   if( !ref $self || ! $self->isa('Bio::EnsEMBL::DB::VirtualContigI') ) {
       $self->throw("Can only extend a VirtualContigI, Bailing out...");
   }
   # based on an original idea by Ewan Birney. ;)
   my $nvc = $self->extend(0,10000000000);
   return $nvc;
}

=head2 windowed_VirtualContig

 Title   : windowed_VirtualContig
 Usage   : $newvc = $vc->windowed_VirtualContig(13400,0,2000);
 Function: Provides a new vc from a current vc as a window
           on the old vc. 
 Example :
 Returns : 
 Args    : position in vc, right distance, left distance.


=cut

sub windowed_VirtualContig {
    my ($self,$position,$left,$right) = @_;

    if( $position < 0 || $position > $self->length ) {
       $self->throw("Attempting to build a new virtual contig out of length bounds!");
    }

    # scan along right->left until we find the first contig
    # whoes start point is before the position
    # hmm - should we presort these guys sometime?
    # Tony: maybe - this module does 3 sorts on the same array...

    my @ids = keys %{$self->{'contighash'}};
    @ids = sort { $self->{'start'}->{$b} <=> $self->{'start'}->{$a} } @ids;
    print STDERR "In Windowed vc...\n";

    print STDERR @ids."\n\n";
    
    my $id = undef;
    foreach $_ ( @ids ) {
	print STDERR "Looking at $_\n";

       if( $self->start_in_vc($_) < $position ) {
           print STDERR "Starting contig: $_ [". $self->start_in_vc($_). " < $position] \n";
           $id = $_;
           last;
       }
    }

    # $id is going to be our new focus. Now - just call 
    # a constructor with appropriate arithmetic...
    my $rc  = $self->{'contighash'}->{$id};
    my $ori = $self->ori_in_vc($id);

    my $wvcpos;
    if( $ori == 1 ) {
       $wvcpos = $rc->golden_start + ($position - $self->start_in_vc($id));
    } else {
       $wvcpos = $rc->golden_end   - ($position - $self->start_in_vc($id));
    }

    
    return Bio::EnsEMBL::DB::VirtualContig->new(-focuscontig   => $rc,
					        -focusposition => $wvcpos,
					        -ori           => $ori,
					        -left          => $left,
					        -right         => $right
					        );


}




=head2 primary_seq

 Title   : primary_seq
 Usage   : $seq = $contig->primary_seq();
 Function: Gets a Bio::PrimarySeqI object out from the contig
 Example :
 Returns : Bio::PrimarySeqI object
 Args    :


=cut

sub primary_seq {
   my ($self) = @_;

   my $seq = $self->_seq_cache();

   if( defined $seq ) {
       return $seq;
   }

   # we have to move across the map, picking up the sequences,
   # truncating them and then adding them into the final product.

   my @contig_id = sort { $self->{'start'}->{$a} <=> $self->{'start'}->{$b} } keys %{$self->{'start'}};

   my $seq_string;
   my $last_point = 1;

   # if there is a left overhang, add it 

   if( $self->_left_overhang() > 0 ) {
       $seq_string = 'N' x $self->_left_overhang();
   }

   foreach my $cid ( @contig_id ) {
       #print(STDERR "\nFinding sequence for $cid\n");
       my $c    = $self->{'contighash'}->{$cid};
       my $tseq = $c->primary_seq();

       #print(STDERR "Seq length/start is " . $tseq->length . "\t" .$self->{start}{$cid} . "\n");
       if( $self->{'start'}->{$cid} != ($last_point+1) ) {
       
           # Tony: added a throw here - if we get negative numbers of inserted N's
	   # my $no = $self->{'start'}->{$cid} - $last_point -1;

	   my $no = $self->{'start'}->{$cid} - $last_point;

           if ($no < 0){
	       $self->throw("Error. Trying to insert negative number ($no) of N\'s into contig sequence");
           }
	   
	   #print STDERR "Putting in $no x N\n";

	   $seq_string .= 'N' x $no;
	   $last_point += $no;
       } 

       my $trunc;
       my $end;

       if( $self->_clone_map == 1 ) {
	   $end = $c->length;
       } else {
	   if( $self->{'rightmostcontig_id'} eq $cid ) {
	       print STDERR "Rightmost end is ",$self->{'rightmostend'},"\n";
	       
	       $end = $self->{'rightmostend'};
	   } else {
	       if( $self->{'contigori'}->{$cid} == 1 ) {
		   $end = $c->golden_end;
	       } else {
		   $end = $c->golden_start;
	       }
	   }
       }
       
       
       print STDERR "got $cid, from ",$self->{'startincontig'}->{$cid}," to ",$c->golden_end,"\n";
        
       if( $self->{'contigori'}->{$cid} == 1 ) {
	   print(STDERR "ori " . $self->{'contigori'}->{$cid} . "\t" . $self->{'startincontig'}->{$cid} . "\t" . $end . "\n");
	   $trunc = $tseq->subseq($self->{'startincontig'}->{$cid},$end);
       } else {
	   print(STDERR "ori " . $self->{'contigori'}->{$cid} . "\t" . $self->{'startincontig'}->{$cid} . "\t" . $end . "\n");
	   $trunc = $tseq->trunc($end,$self->{'startincontig'}->{$cid})->revcom->seq;
       }
#       print STDERR "Got $trunc\n";
       $seq_string .= $trunc;
       $last_point += length($trunc);
   }

   # if there is a right overhang, add it 

   if( $self->_right_overhang() > 0 ) {
       $seq_string .= 'N' x $self->_right_overhang();
   }


   $seq = Bio::PrimarySeq->new( -id      => "virtual_contig_".$self->_unique_number,
				-seq     => $seq_string,
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

sub id {
    my ($self) = @_;

    return "virtual_contig_".$self->_unique_number;
}

=head2 top_SeqFeatures

 Title   : top_SeqFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub top_SeqFeatures {
   my ($self,@args) = @_;
   my (@f);


   if( $self->skip_SeqFeature('similarity') != 1 ) { 
       push(@f,$self->get_all_SimilarityFeatures());
   } 

   if( $self->skip_SeqFeature('repeat') != 1 ) { 
       push(@f,$self->get_all_RepeatFeatures());
   } 

   if( $self->skip_SeqFeature('external') != 1 ) { 
       push(@f,$self->get_all_ExternalFeatures());
   } 

   foreach my $gene ( $self->get_all_Genes()) {
       print STDERR "Got a $gene\n";
       my $vg = Bio::EnsEMBL::VirtualGene->new(-gene => $gene,-contig => $self);
       push(@f,$vg);
   }

   return @f;
}


=head2 get_all_SeqFeatures

 Title   : get_all_SeqFeatures
 Usage   : foreach my $sf ( $contig->get_all_SeqFeatures ) 
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_SeqFeatures {
   my ($self) = @_;
   my @out;
   push(@out,$self->get_all_SimilarityFeatures());
   push(@out,$self->get_all_RepeatFeatures());

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

sub get_all_SimilarityFeatures {
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

sub get_all_RepeatFeatures {
   my ($self) = @_;
   
   return $self->_get_all_SeqFeatures_type('repeat');

}

=head2 get_all_ExternalFeatures

 Title   : get_all_ExternalFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_ExternalFeatures {
   my ($self) = @_;

   return $self->_get_all_SeqFeatures_type('external');
}

=head2 get_all_PredictionFeatures

 Title   : get_all_PredictionFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_PredictionFeatures {
   my ($self) = @_;

   return $self->_get_all_SeqFeatures_type('prediction');
}



=head2 get_all_Genes

 Title   : get_all_Genes
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_Genes {
    my ($self) = @_;
    my (%gene,%trans,%exon);

    foreach my $c ( values %{$self->{'contighash'}} ) {
	foreach my $gene ( $c->get_all_Genes() ) {
	    $gene{$gene->id()} = $gene;
	}
    }
    
    # get out unique set of translation objects
    foreach my $gene ( values %gene ) {
	foreach my $transcript ( $gene->each_Transcript ) {
	    my $translation = $transcript->translation;
	    $trans{"$translation"} = $translation;
	    
	}
    }

    foreach my $gene ( values %gene ) {
	foreach my $exon ( $gene->all_Exon_objects() ) {
	    # hack to get things to behave
	    #print STDERR "Exon was ",$exon->start,":",$exon->end,":",$exon->strand,"\n";
	    $exon->seqname($exon->contig_id);
	    $exon{$exon->id} = $exon;
	    $self->_convert_seqfeature_to_vc_coords($exon);
	    #print STDERR "Exon going to ",$exon->start,":",$exon->end,":",$exon->strand," ,",$exon->seqname,"\n";
	}
    }
    
    foreach my $t ( values %trans ) {
	if( exists $self->{'contighash'}->{$exon{$t->start_exon_id}->contig_id} ) {
	    my ($start,$end,$str) = $self->_convert_start_end_strand_vc($exon{$t->start_exon_id}->contig_id,$t->start,$t->start,1);
	    $t->start($start);
	}

	if( exists $self->{'contighash'}->{$exon{$t->end_exon_id}->contig_id} ) {
	    my ($start,$end,$str) = $self->_convert_start_end_strand_vc($exon{$t->end_exon_id}->contig_id,$t->end,$t->end,1);
	    $t->end($start);
	}
    }
    
    return values %gene;
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

   #return $self->_left_size + $self->_right_size +1;
   return $self->_left_size + $self->_right_size;

}

=head2 embl_accession

 Title   : embl_accession
 Usage   : $obj->embl_accession($newval)
 Function: 
 Returns : value of embl_accession
 Args    : newvalue (optional)


=cut

sub embl_accession {
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'embl_accession'} = $value;
    }
    return $obj->{'embl_accession'};

}


=head2 dbobj

 Title   : dbobj
 Usage   : $obj->dbobj($newval)
 Function: 
 Returns : value of dbobj
 Args    : newvalue (optional)


=cut

sub dbobj {
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'dbobj'} = $value;
    }
    return $obj->{'dbobj'};

}

=head2 skip_SeqFeature

 Title   : skip_SeqFeature
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub skip_SeqFeature {
   my ($self,$tag,$value) = @_;

   if( defined $value ) {
       $self->{'feature_skip'}->{$tag} = $value;
   }

   return $self->{'feature_skip'}->{$tag};
}

=head2 rawcontig_ids

 Title   : rawcontig_ids
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub rawcontig_ids {
   my ($self,@args) = @_;

   return keys %{$self->{'contighash'}};
}

=head2 start_in_vc

 Title   : start_in_vc
 Usage   : $vc->start_in_vc('rawcontigid');
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub start_in_vc {
   my ($self,$rawcontigid) = @_;
   
   if( !defined $rawcontigid || ! exists $self->{'contighash'}->{$rawcontigid} ) {
       $self->throw("No rawcontig id provided/not in vc [$rawcontigid]");
   }

   return $self->{'start'}->{$rawcontigid}
}

=head2 end_in_vc

 Title   : end_in_vc
 Usage   : $vc->end_in_vc('rawcontigid')
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub end_in_vc {
   my ($self,$cid) = @_;

   if( $self->{'rightmostcontigid'} eq $cid ) {
       return $self->start_in_vc($cid) + $self->{'rightmostend'};
   } else {
       return $self->start_in_vc($cid) + $self->{'contighash'}->{$cid}->golden_length;
   } 

}

sub ori_in_vc {
    my ($self,$rawcontigid) = @_;
    
    if( !defined $rawcontigid || ! exists $self->{'contighash'}->{$rawcontigid} ) {
	$self->throw("No rawcontig id provided/not in vc [$rawcontigid]");
    }
    
    return $self->{'contigori'}->{$rawcontigid}
}



=head2 _build_clone_map

 Title   : _build_clone_map
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _build_clone_map {
   my ($self,$clone) = @_;

   my $total_len   = 0;
   my $length      = 0;
   my $seen        = 0;
   my $middle      = 0;
   
   #print STDERR "Making clone map\n";
   foreach my $contig ( $clone->get_all_Contigs ) {
       $self->{'start'}        ->{$contig->id} = $contig->embl_offset;
       $self->{'startincontig'}->{$contig->id} = 1;
       $self->{'contigori'}    ->{$contig->id} = 1;
       $self->{'contighash'}   ->{$contig->id} = $contig;

       #print STDERR "Clone map: contig: ",$contig->id," [em_offset: ",$contig->embl_offset,"] ==> [contig_offset: ",$contig->length,"]\n";

       $total_len = $contig->embl_offset + $contig->length;

       if( $total_len > $length ) {
	   $length = $total_len;
           #print STDERR "New contig length: ",$length,"\n";
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
   $self->_left_size($middle-1);
   $self->_right_size($length-$middle);

   # Remember this vc contructed from a clone (rather than extending a 'seed' contig)
   $self->_clone_map(1);

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

sub _build_contig_map {
  my ($self,$focuscontig,$focusposition,$ori,$left,$right) = @_;
  
  # we first need to walk down contigs going left
  # so we can figure out the start position (contig-wise)
  
  # initialisation - find the correct end of the focus contig
  print(STDERR "in build_contig_map\n");
  my ($current_left_size,$current_orientation,$current_contig,$overlap);
    
  if( $ori == 1 ) {
    $current_left_size   = $focusposition;
    $current_orientation = 1;
  } else {
    $current_left_size   = $focuscontig->length - $focusposition;
    $current_orientation = -1;
  }
  
  $current_contig = $focuscontig;
  
  print STDERR "Left size is $left\n";
  
  GOING_LEFT :
    
    while( $current_left_size < $left ) {
      print(STDERR "Current left = $current_left_size\n");
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
	
	if( $overlap->distance == 1 ) {
	  $current_left_size += $overlap->sister->golden_length -1;
	} else {
	  $current_left_size += $overlap->distance;
	  if( $current_left_size > $left ) {
	    # set the left overhang!
	    print STDERR "Triggered left overhang - ",$overlap->distance,":$current_left_size\n";
	    $self->_left_overhang($overlap->distance - ($current_left_size - $left));
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
	  last;
	}
	
	if( $overlap->distance == 1 ) {
	  $current_left_size += $overlap->sister->golden_length-1;
	  print STDERR ("Current " . $current_left_size . " " . $left . "\n");
	} else {
	  $current_left_size += $overlap->distance;
	  print STDERR ("Current " . $current_left_size . " " . $left . "\n");
	  if( $current_left_size > $left ) {
	    # set the left overhang!
	    print(STDERR "Setting left overhang " . $current_left_size . " " . $left . "\n");
	    $self->_left_overhang($overlap->distance - ($current_left_size - $left));
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
    }
  
  # now $current_contig is the left most contig in this set, with
  # its orientation set and ready to rock... ;)
  
  my $total = $left + $right;
  
  print STDERR "leftmost contig is "  . $current_contig->id . 
               " with $total to account for, " .
	       " gone $current_left_size of $left\n";

  $self->{'leftmostcontig_id'} = $current_contig->id;
  
  # the first contig will need to be trimmed at a certain point
  my $startpos;
  
  print STDERR "Leftmost contig starts at: $startpos orientation: $current_orientation\n";
  print STDERR "Current left = $left vs global left= $current_left_size\n";
  
  my $current_length;
  print STDERR "Left overhang is " . $self->_left_overhang . "\n";
  if( $self->_left_overhang() == 0 ) {
    if( $current_orientation == 1 ) {
      print(STDERR "Current orientation $current_orientation\n");
      print(STDERR "golden start " . $current_contig->golden_start . "\n");
      $startpos = $current_contig->golden_start + ($current_left_size - $left);
    } else {
      print(STDERR "Current orientation $current_orientation\n");
      print(STDERR "golden end " . $current_contig->golden_end . "\n");
      
      $startpos = $current_contig->golden_end   - ($current_left_size - $left);
    }
	
	print STDERR "Leftmost contig has $startpos and $current_orientation $left vs $current_left_size\n";
	
	$self->{'start'}        ->{$current_contig->id} = 1;
	$self->{'startincontig'}->{$current_contig->id} = $startpos;
	$self->{'contigori'}    ->{$current_contig->id} = $current_orientation;
	$self->{'contighash'}   ->{$current_contig->id} = $current_contig;

	if( $current_orientation == 1 ) {
	  $current_length = $current_contig->golden_end - $startpos +1;
	} else {
	  $current_length = $startpos - $current_contig->golden_start+1;
	}
  } else {
	# has an overhang - first contig offset into the system
    print STDERR "First contig offset " . $self->_left_overhang . "\n";
      $self->{'start'}        ->{$current_contig->id} = $self->_left_overhang+1;
	if( $current_orientation == 1 ) {
	  $self->{'startincontig'}->{$current_contig->id} = $current_contig->golden_start;
	} else {
	  print STDERR "Start in contig is " . $current_contig->golden_end . "\n";

	  $self->{'startincontig'}->{$current_contig->id} = $current_contig->golden_end;
	}
	$self->{'contigori'}    ->{$current_contig->id} = $current_orientation;
	$self->{'contighash'}   ->{$current_contig->id} = $current_contig;

	$current_length = $self->_left_overhang() + $current_contig->golden_length ;
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
		$self->found_right_end(1);
		$right = $current_length - $left;
		last;
	   }
	    
	    # see whether the distance gives us an end condition, and a right_overhang

	    if( $current_length + $overlap->distance > $total ) {
		# right overhang
		$self->_right_overhang($total - $current_length);
		last;
	    }

	    # add to total, move on the contigs
	    
	    $current_contig = $overlap->sister();
	    print(STDERR "New sister " . $current_contig->id  . "\n");
	    $self->{'contighash'}->{$current_contig->id} = $current_contig;
	    
	    if( $overlap->sister_polarity == 1) {
		$current_orientation = 1;
	    } else {
		$current_orientation = -1;
	    }
	    
	    # The +1's here are to handle the fact we want to produce abutting
	    # coordinate systems from overlapping switch points.
	    if( $overlap->distance == 1 ) {
		$self->{'start'}->{$current_contig->id} = $current_length +1;
	    } else {
		$self->{'start'}->{$current_contig->id} = $current_length + $overlap->distance;
		$current_length += $overlap->distance;
	    }
	    
	    if( $current_orientation == 1 ) {
		if( $overlap->distance == 1 ) {
		    $self->{'startincontig'}->{$current_contig->id} = $current_contig->golden_start+1; 
		} else {
		    $self->{'startincontig'}->{$current_contig->id} = $current_contig->golden_start;
		}
	    } else {
		if( $overlap->distance == 1 ) {
		    $self->{'startincontig'}->{$current_contig->id} = $current_contig->golden_end-1; 
		} else {
		    $self->{'startincontig'}->{$current_contig->id} = $current_contig->golden_end; 
		}
	    }
	    
	    $self->{'contigori'}->{$current_contig->id} = $current_orientation;
	    
	    # up the length
	    $current_length += $overlap->sister->golden_length -1;
	} else {
	    # go left wrt to the contig
	    print STDERR "Going left with " . $current_contig->id . "\n";
	    
	    $overlap = $current_contig->get_left_overlap();
	 
	    # if there is no left overlap, trim right to this size
	    # as this means we have run out of contigs
	    if( !defined $overlap ) {
		$self->found_left_end(1);
		$right = $current_length - $left;
		last;
	    }

	   # add to total, move on the contigs

	   $current_contig = $overlap->sister();
	    print(STDERR "New sister " . $current_contig->id  . "\n");
	   $self->{'contighash'}->{$current_contig->id} = $current_contig;

	   if( $overlap->sister_polarity == 1) {
	       $current_orientation = -1;
	   } else {
	       $current_orientation = 1;
	   }

	   # The +1's here are to handle the fact we want to produce abutting
	   # coordinate systems from overlapping switch points.

	   if( $overlap->distance == 1 ) {
	       $self->{'start'}->{$current_contig->id} = $current_length +1;
	   } else {
	       $self->{'start'}->{$current_contig->id} = $current_length + $overlap->distance;
	       $current_length += $overlap->distance;
	   }

	   if( $current_orientation == 1 ) {
	       $self->{'startincontig'}->{$current_contig->id} = $current_contig->golden_start +1;
	   } else {
	       $self->{'startincontig'}->{$current_contig->id} = $current_contig->golden_end -1;
	   }

	   $self->{'contigori'}->{$current_contig->id} = $current_orientation;

	   # up the length
	   $current_length += $overlap->sister->golden_length -1;
       }
   }

   # $right might have been modified during the walk

   $total = $left + $right;

   # need to store end point for last contig
   print STDERR "Looking at setting rightmost end with $total and $current_length ",$current_contig->golden_end,"\n";

    $self->{'rightmostcontig_id'} = $current_contig->id();
    if( $self->_right_overhang == 0 ) {
	if( $current_orientation == 1 ) {
	    $self->{'rightmostend'}    = $current_contig->golden_end - ($current_length - $total);
	} else {
	    $self->{'rightmostend'}    = $current_contig->golden_start + ($current_length - $total);
	}
    } else {
	if( $current_orientation == 1 ) {
	    $self->{'rightmostend'}    = $current_contig->golden_end;
	} else {
	    $self->{'rightmostend'}    = $current_contig->golden_start;
	}
    }
   
   # put away the focus/size info etc

   $self->_focus_contig($focuscontig);
   $self->_focus_position($focusposition);
   $self->_focus_orientation($ori);
   $self->_left_size($left);
   $self->_right_size($right);


   # ready to rock and roll. Woo-Hoo!

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

	$self->{_found_left_end} = $arg;
    } elsif (defined($arg)) {
	$self->throw("Arg to found_left_end should be 0,1");
    }

    return $self->{_found_left_end};
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
	$self->{_found_right_end} = $arg;
    } elsif (defined($arg)) {
	$self->throw("Arg to found_right_end should be 0,1");
    }

    return $self->{_found_right_end};
}

=head2 is_truncated

 Title   : is_truncated
 Usage   : 
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub is_truncated {
    my ($self) = @_;

    my $flag = 0;

    if ($self->found_right_end == 1) {
	$flag = 1;
    } 
    if ($self->found_left_end == 1) {
	$flag = 1;
    }

    return $flag;
}

=head2 _get_all_SeqFeatures_type

 Title   : _get_all_SeqFeatures_type
 Usage   : Internal function which encapsulates getting
           features of a particular type and returning
           them in the VC coordinates, optionally cacheing
           them.
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _get_all_SeqFeatures_type {
   my ($self,$type) = @_;

   if( $self->_cache_seqfeatures() && $self->_has_cached_type($type) ) {
       return $self->_get_cache($type);
   }

   # ok - build the sequence feature list...

   my $sf;

   if( $self->_cache_seqfeatures() ) {
       $sf = $self->_make_cache($type);
   } else {
       $sf = []; 
   }

   foreach my $c ( values %{$self->{'contighash'}} ) {
       #print STDERR "Looking at ",$c->id,"\n";
       if( $type eq 'repeat' ) {
	   push(@$sf,$c->get_all_RepeatFeatures());
       } elsif ( $type eq 'similarity' ) {
	   push(@$sf,$c->get_all_SimilarityFeatures());
       } elsif ( $type eq 'prediction' ) {
	   push(@$sf,$c->get_all_PredictionFeatures());
       } elsif ( $type eq 'external' ) {
	   push(@$sf,$c->get_all_ExternalFeatures());
       } else {
	   $self->throw("Type $type not recognised");
       }
   }

   my @vcsf = ();
   # need to clip seq features to fit the boundaries of
   # our v/c so displays don't break

   my $count = 0;
   foreach $sf ( @$sf ) {
       $sf = $self->_convert_seqfeature_to_vc_coords($sf);
       if( !defined $sf ) {
	   next;
       }

	
       if($sf->start < 0 ){
            #print STDERR "Discarding (L) vc feature: ",$sf->seqname," at start: ",$sf->start," end: ",$sf->end,"\n";
            $count++;
        }
        elsif ($sf->end > $self->length){
            #print STDERR "Discarding (R) vc feature: ",$sf->seqname," at start: ",$sf->start," end: ",$sf->end,"\n";
            $count++;
        }
        else{
            #print STDERR "Keeping vc feature: ",$sf->seqname," at start: ",$sf->start," end: ",$sf->end,"\n";
            push (@vcsf, $sf);
        }
   }
   
   return @vcsf;
}


=head2 _convert_seqfeature_to_vc_coords

 Title   : _convert_seqfeature_to_vc_coords
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _convert_seqfeature_to_vc_coords {
   my ($self,$sf) = @_;

   my $cid = $sf->seqname();

   if( !defined $cid ) {
       $self->throw("sequence feature [$sf] has no seqname!");
   }

   if( !exists $self->{'contighash'}->{$cid} ) {
       return undef;
   }

   # if this is something with subfeatures, then this is much more complex
   my @sub = $sf->sub_SeqFeatures();

   if( $#sub >=  0 ) {
       # chain to constructor of the object. Not pretty this.
       my $new = $sf->new();

       if( $new->can('attach_seq') ) {
	   $new->attach_seq($self->primary_seq);
       }

       my $seen = 0;
       foreach my $sub ( @sub ) {
	   $sub = $self->_convert_seqfeature_to_vc_coords($sub);
	   if( !defined $sub ) {
	       next;
	   }
	   $seen =1;
	   $new->add_SeqFeature($sub);
       }
       if( $seen == 1 ) {
	   return $new;
       } else {
	   return undef;
       }
   }



   # might be clipped left/right

   if( $self->{'leftmostcontig_id'} eq $cid ){
       if( $self->ori_in_vc($cid) == 1) {
	   # if end is less than startincontig - a no-go
	   if( $sf->end < $self->{'startincontig'}->{$cid} ) {
	       return 0;
	   }
       } else {
	   # if start is > start in contig
	   if( $sf->start > $self->{'startincontig'}->{$cid} ) {
	       return 0;
	   }
       }
   }


   my ($rstart,$rend,$rstrand) = $self->_convert_start_end_strand_vc($cid,$sf->start,$sf->end,$sf->strand);
   
   $sf->start ($rstart);
   $sf->end   ($rend);
   $sf->strand($rstrand);

   if( $sf->can('attach_seq') ) {
       $sf->attach_seq($self->primary_seq);
   }

   $sf->seqname($self->id);

   return 1;
}

=head2 _convert_start_end_strand_vc

 Title   : _convert_start_end_strand_vc
 Usage   : Essentially an internal for _convert_seqfeature,
           but sometimes we have coordiates  not on seq features
 Function:
 Example : ($start,$end,$strand) = $self->_convert_start_end_strand_vc($contigid,$start,$end,$strand)
 Returns : 
 Args    :


=cut

sub _convert_start_end_strand_vc {
   my ($self,$contig,$start,$end,$strand) = @_;
   my ($rstart,$rend,$rstrand);

   if( !exists $self->{'contighash'}->{$contig} ) {
       $self->throw("Attempting to map a sequence feature with [$contig] on a virtual contig with no $contig");
   }
   if( $self->{'contigori'}->{$contig} == 1 ) {
       # ok - forward with respect to vc. Only need to add offset
       my $offset = $self->{'start'}->{$contig} - $self->{'startincontig'}->{$contig};
       $rstart = $start + $offset;
       $rend   = $end + $offset;
       $rstrand = $strand;
   } else {
       my $offset = $self->{'start'}->{$contig} + $self->{'startincontig'}->{$contig};
       # flip strand
       $rstrand = $strand * -1;
       
       # yup. A number of different off-by-one errors possible here


       $rstart  = $offset - $end;
       $rend    = $offset - $start;

   }

   return ($rstart,$rend,$rstrand);
}



=head2 _dump_map

 Title   : _dump_map
 Usage   : Produces a dumped map for debugging purposes
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _dump_map {
   my ($self,$fh) = @_;

   ! defined $fh && do { $fh = \*STDERR};
 
   my @ids = keys %{$self->{'contighash'}};
   @ids = sort { $self->{'start'}->{$a} <=> $self->{'start'}->{$b} } @ids;

   print $fh "Contig Map Dump: \n";
   foreach my $id ( @ids ) {
       print $fh "Contig $id starts:",$self->{'start'}->{$id}," start in contig ",$self->{'startincontig'}->{$id}," orientation ",$self->{'contigori'}->{$id},"\n";
   }
}


=head2 _focus_contig

 Title   : _focus_contig
 Usage   : $obj->_focus_contig($newval)
 Function: 
 Example : 
 Returns : value of _focus_contig
 Args    : newvalue (optional)


=cut

sub _focus_contig {
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'_focus_contig'} = $value;
    }
    return $obj->{'_focus_contig'};

}

=head2 _focus_position

 Title   : _focus_position
 Usage   : $obj->_focus_position($newval)
 Function: 
 Example : 
 Returns : value of _focus_position
 Args    : newvalue (optional)


=cut

sub _focus_position {
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

sub _focus_orientation {
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

sub _left_size {
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

sub _right_size {
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

sub _cache_seqfeatures {
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

sub _has_cached_type {
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

sub _make_cache {
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

sub _get_cache {
   my ($self,$type) = @_;

   return $self->{'_sf_cache'}->{$type};
   
}

=head2 _clear_vc_cache

 Title   : _clear_vc_cache
 Usage   : $virtual_contig->_clear_vc_cache
 Function: clears a v/c internal cache (use when extending a virtual contig)
 Example :
 Returns : 
 Args    :


=cut

sub _clear_vc_cache {
    my $self = shift;
    $self->{'_seq_cache'} = undef;
    foreach my $c (keys %{$self->{'_sf_cache'}}){
        $self->{'_sf_cache'}->{$c} = undef;
    }   
}

=head2 _unique_number

 Title   : _unique_number
 Usage   : $obj->_unique_number($newval)
 Function: 
 Returns : value of _unique_number
 Args    : newvalue (optional)


=cut

sub _unique_number {
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

sub _seq_cache {
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_seq_cache'} = $value;
    }
    return $obj->{'_seq_cache'};

}

=head2 _clone_map

 Title   : _clone_map
 Usage   : $obj->_clone_map($newval)
 Function: 
 Example : 
 Returns : value of _clone_map
 Args    : newvalue (optional)


=cut

sub _clone_map {
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'_clone_map'} = $value;
    }
    return $obj->{'_clone_map'};

}


=head2 _at_left_end

 Title   : _at_left_end
 Usage   : $obj->_at_left_end($newval)
 Function: 
 Example : 
 Returns : true if we have reached the  obsolute left end of vc
 Args    : newvalue (optional)


=cut


sub _at_left_end {
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'_at_left_end'} = $value;
    }
    return $obj->{'_at_left_end'};

}

=head2 _at_right_end

 Title   : _at_right_end
 Usage   : $obj->_at_right_end($newval)
 Function: 
 Example : 
 Returns : true if we have reached the  obsolute right end of vc
 Args    : newvalue (optional)


=cut


sub _at_right_end {
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'_at_right_end'} = $value;
    }
    return $obj->{'_at_right_end'};

}

=head2 _left_overhang

 Title   : _left_overhang
 Usage   : $obj->_left_overhang($newval)
 Function: 
 Example : 
 Returns : value of _left_overhang
 Args    : newvalue (optional)


=cut

sub _left_overhang{
   my ($obj,$value) = @_;
   
   if( defined $value) {
     print(STDERR "Setting overhand to $value\n");
      $obj->{'_left_overhang'} = $value;
    }
    return $obj->{'_left_overhang'};

}

=head2 _right_overhang

 Title   : _right_overhang
 Usage   : $obj->_right_overhang($newval)
 Function: 
 Example : 
 Returns : value of _right_overhang
 Args    : newvalue (optional)


=cut

sub _right_overhang{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'_right_overhang'} = $value;
    }
    return $obj->{'_right_overhang'};

}


1;






