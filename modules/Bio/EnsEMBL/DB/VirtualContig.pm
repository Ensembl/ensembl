
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
use Bio::EnsEMBL::DB::VirtualMap;
use Bio::EnsEMBL::DB::VirtualPrimarySeq;

my $VC_UNIQUE_NUMBER = 0;

@ISA = qw(Bio::Root::Object Bio::EnsEMBL::DB::VirtualContigI);

# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;
  
  my $make = $self->SUPER::_initialize(@args);
  
  $self->name("Virtual Contig Module"); # set the exception context (does this work?)

  my ($focuscontig,$focusposition,$ori,$leftsize,$rightsize,$clone) = $self->_rearrange([qw( FOCUSCONTIG FOCUSPOSITION ORI LEFT RIGHT CLONE)],@args);
  
  #Create a new VirtualMap holder object for MapContigs
  my $vmap=Bio::EnsEMBL::DB::VirtualMap->new();
  $self->_vmap($vmap);

  # set up feature_skip hash
  $self->{'feature_skip'}  = {};

  # this is for cache's of sequence features if/when we want them
  $self->{'_sf_cache'}     = {};

  $self->_vmap->left_overhang(0);
  $self->_vmap->right_overhang(0);

  if( defined $clone && defined $focuscontig ){
      $self->throw("Build a virtual contig either with a clone or a focuscontig, but not with both");
  }
  
  if( defined $clone ) {
      $self->_vmap->build_clone_map($clone);
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
      if (! $focuscontig->isa('Bio::EnsEMBL::DBSQL::RawContig') ) {
	  $self->throw("$focuscontig is not a Bio::EnsEMBL::DBSQL::RawContig object, cannot make Virtual Contig!");
      }
      
      # build the map of how contigs go onto the vc coorindates
      $self->_vmap->build_contig_map($focuscontig,$focusposition,$ori,$leftsize,$rightsize);
      $self->_vmap->dbobj($focuscontig->dbobj);
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
    
    # checky $self call to make sure we get the original type of object...

    my $nvc = $self->new( -focuscontig => $self->_vmap->focus_contig,
						    -focusposition   => $self->_vmap->focus_position,
						    -ori             => $self->_vmap->focus_orientation,
						    -left            => $self->_vmap->left_size - $left,
						    -right           => $self->_vmap->right_size + $right,
						    );
    
    my $id = join('.', ($nvc->_vmap->focus_contig->id, $nvc->_vmap->focus_position, $nvc->_vmap->focus_orientation, $nvc->_vmap->left_size, $nvc->_vmap->right_size));
    $nvc->_unique_number($id);
    
    return $nvc;
}


sub get_MarkerFeatures {

    my ($self)=@_;
     return $self->_get_all_SeqFeatures_type('marker');
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
    my $nvc = $self->extend(-10000000000,10000000000);
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
 Args    : focus position in vc For wvc, left distance to extend, right distance to extend.


=cut

sub windowed_VirtualContig {
    my ($self,$position,$left,$right) = @_;
    
    if( $position < 0 || $position > $self->length ) {
	$self->throw("Attempting to build a new virtual contig out of length bounds!");
    }
    
    my ($map_contig,$f_position,$ori) = $self->raw_contig_position($position,1);
    
    return Bio::EnsEMBL::DB::VirtualContig->new(-focuscontig   => $map_contig,
					        -focusposition => $f_position,
					        -ori           => $ori,
					        -left          => $left,
					        -right         => $right
					        );
}

=head2 virtual_primary_seq

 Title   : virtual_primary_seq
 Usage   : $seq = $contig->virtual_primary_seq();
 Function: Gets a Bio::EnsEMBL::DB::VirtualPrimarySeq object out from the contig
 Returns : Bio::EnsEMBL::VirtualPrimarySeq object
 Args    : none

=cut

sub virtual_primary_seq {
    my ($self) = @_;
    
    my $vseq = Bio::EnsEMBL::DB::VirtualPrimarySeq->new(-vmap =>$self->_vmap,
						 -un =>$self->_unique_number,
						 -clone =>$self->_vmap->clone_map,
						 -"length" =>$self->length
						 );
    return $vseq;
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

   my @map_contigs=$self->_vmap->get_all_MapContigs;

   my $seq_string;
   my $last_point = 1;
   
   # if there is a left overhang, add it 
   
   if( $self->_vmap->left_overhang() > 0 ) {
       $seq_string = 'N' x $self->_vmap->left_overhang();
   }
   
   #Go through each MapContig
   my $previous = undef;
   foreach my $mc ( @map_contigs ) {
       if( defined $previous && $previous->end+1 != $mc->start ) {
	   # then start had better be before end
	   if( $mc->start < $previous->end ) {
	       $self->throw("Inconsistent map contigs with start of next contig less than end of previous");
	   }

	   my $length = $mc->start - $previous->end -1;
	   my $str = 'N' x $length;
	   #print STDERR "Adding in $length N's\n";
	   $seq_string .= $str;
       } else {
	   # exact overlap - chew back the switch base
	   #$seq_string = substr($seq_string,0,-1);
	   # chew back no longer required...
       }
       
       # now add in the actual sequence.
       #print STDERR $mc->start," Start in is ".$mc->contig->id." ".$mc->start_in.":".$mc->end_in." ".$mc->contig->length."\n";
       #print STDERR "adding to",length($seq_string),"\n";

       # flip if other way around
       my $substring;
       if( $mc->orientation == -1 ) {
	   $substring = $mc->contig->primary_seq->subseq($mc->end_in,$mc->start_in);
	   $substring =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
	   $substring = CORE::reverse $substring;
       } else {
	   $substring = $mc->contig->primary_seq->subseq($mc->start_in,$mc->end_in);
       }
       
       # add it
       $seq_string .= $substring;
       $previous = $mc;
   }
   
   # if there is a right overhang, add it 
   if( $self->_vmap->right_overhang() > 0 ) {
       $seq_string .= 'N' x $self->_vmap->right_overhang();
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

    
    if( !$self->skip_SeqFeature('similarity')  ) { 
	push(@f,$self->get_all_SimilarityFeatures());
    } 
    
    if( !$self->skip_SeqFeature('repeat')  ) { 
	push(@f,$self->get_all_RepeatFeatures());
    } 
    
    if( !$self->skip_SeqFeature('external')  ) { 
	push(@f,$self->get_all_ExternalFeatures());
    } 
    
    foreach my $gene ( $self->get_all_Genes()) {
	my $vg = Bio::EnsEMBL::VirtualGene->new(-gene => $gene,-contig => $self);
	push(@f,$vg);
    }
    
    return @f;
}

=head2 get_all_VirtualGenes
    
 Title   : get_all_VirtualGenes
 Usage   : foreach my $virtualgene ( $contig->get_all_VirtualGenes ) 
 Function: Gets all teh genes on this VC as VirtualGene objects
 Example : 
 Returns : array of Bio::EnsEMBL::VirtualGene objects
 Args    : none


=cut

sub get_all_VirtualGenes {
    my ($self) = @_;
    
    my @out;
    foreach my $gene ( $self->get_all_Genes()) {
	my $vg = Bio::EnsEMBL::VirtualGene->new(-gene => $gene,-contig => $self);
	push(@out,$vg);
    }
    return @out;
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


=head2 get_all_SimilarityFeatures_above_score

 Title   : get_all_SimilarityFeatures_above_score
 Usage   : foreach my $sf ( $contig->get_all_SimilarityFeatures_above_score(analysis_type, score) ) 
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_SimilarityFeatures_above_score{
    my ($self, $analysis_type, $score) = @_;
    
    my $sf = ();
    foreach my $c ($self->_vmap->get_all_RawContigs) {
	   push(@$sf, $c->get_all_SimilarityFeatures_above_score($analysis_type, $score));
   }

   # Need to clip seq features to fit the boundaries of
   # our v/c so displays don't break
   my @vcsf = ();
   my $count = 0;
   foreach $sf ( @$sf ) {

       $sf = $self->_convert_seqfeature_to_vc_coords($sf);

       if( !defined $sf ) {      
	   next;
       }	
       if (($sf->start < 0 ) || ($sf->end > $self->length)) {
	   $count++;
       }
       else{
	    push (@vcsf, $sf);
       }
   }
   
   return @vcsf;
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
    my (%gene,%trans,%exon,%exonconverted);
    
    foreach my $contig ($self->_vmap->get_all_RawContigs) {
	foreach my $gene ( $contig->get_all_Genes() ) {      
	    $gene{$gene->id()} = $gene;
	}
    }

    return $self->_gene_query(%gene);

}




=head2 get_Genes_by_Type

 Title   : get_Genes_by_Type
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Genes_by_Type {
    my ($self,$type) = @_;
    my (%gene,%trans,%exon,%exonconverted);
    
    foreach my $contig ($self->_vmap->get_all_RawContigs) {
	foreach my $gene ( $contig->get_Genes_by_Type($type) ) {      
	    $gene{$gene->id()} = $gene;
	}
    }

    return $self->_gene_query(%gene);

}




sub _gene_query{
    
    my ($self,%gene) = @_;
    my (%trans,%exon,%exonconverted);
        

    foreach my $gene ( values %gene ) {
	
        my $internalExon = 0;
	foreach my $exon ( $gene->all_Exon_objects() ) {
	    # hack to get things to behave
	    $exon->seqname($exon->contig_id);
	    $exon{$exon->id} = $exon;
	    if ($self->_convert_seqfeature_to_vc_coords($exon)) {
                $internalExon = 1;
		$exonconverted{$exon->id} = 1;
            }                           
	}
        
        unless ($internalExon) {    
            delete $gene{$gene->id};
        } 
    }
    
    # get out unique set of translation objects
    foreach my $gene ( values %gene ) {
	foreach my $transcript ( $gene->each_Transcript ) {
	    my $translation = $transcript->translation;
	    $trans{"$translation"} = $translation;	    
	}
    } 
    
    foreach my $t ( values %trans ) {

	if( exists $exonconverted{$t->start_exon_id} ) {
	    my ($start,$end,$str) = $self->_convert_start_end_strand_vc($exon{$t->start_exon_id}->contig_id,$t->start,$t->start,1);
	    $t->start($start);
	}

	if( exists $exonconverted{$t->end_exon_id}  ) {
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

   return $self->_vmap->left_size + $self->_vmap->right_size;
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
 Usage   : my ($map_contig,$rc_position,$rc_strand) = $vmap->raw_contig_position ($vc_pos,$vc_strand)
 Function: Maps a VirtualContig position to the RawContig Position
 Returns : Bio::EnsEMBL::DB::MapContig object, 
           position (int), strand (int)
 Args    : position (int), strand (int)


=cut

sub raw_contig_position {
    my ($self, $vcpos, $vcstrand)=@_;

    return $self->_vmap->raw_contig_position($vcpos,$vcstrand);
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
 Usage   : 
 Function: 
 Returns : value of dbobj
 Args    : 


=cut

sub dbobj {
   my $self =shift;
   return $self->_vmap->dbobj;
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
   my ($self) = @_;
   
   return $self->_vmap->RawContig_ids;
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
    my ($self) = @_;
    return $self->_vmap->found_left_end;
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
    my ($self) = @_;
    
    return $self->_vmap->found_right_end;
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
   
   foreach my $c ($self->_vmap->get_all_RawContigs) {
       if( $type eq 'repeat' ) {
	   push(@$sf,$c->get_all_RepeatFeatures());
       } elsif ( $type eq 'similarity' ) {
	   push(@$sf,$c->get_all_SimilarityFeatures());
       } elsif ( $type eq 'prediction' ) {
	   push(@$sf,$c->get_all_PredictionFeatures());
       } elsif ( $type eq 'external' ) {
	   push(@$sf,$c->get_all_ExternalFeatures());
       } elsif ( $type eq 'marker' ) {
	   push(@$sf,$c->get_MarkerFeatures());
       } else {
	   $self->throw("Type $type not recognised");
       }
   }

   my @vcsf = ();
   # need to clip seq features to fit the boundaries of
   # our v/c so displays don't break

   my $count = 0;
   foreach $sf ( @$sf ) {
       #print "\n ##### Starting to convert featre " . $sf->seqname . " " . $sf->id . "\n";
       $sf = $self->_convert_seqfeature_to_vc_coords($sf);

       if( !defined $sf ) {      
	   next;
       }

	
       if($sf->start < 0 ){
	   $count++;
        }
        elsif ($sf->end > $self->length){
	    $count++;
        }
        else{
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
    my $mc;
    if( !defined $cid ) {
	$self->throw("sequence feature [$sf] has no seqname!");
    }
    eval {
	$mc=$self->_vmap->get_MapContig($cid);
    };
    if ($@) { #print STDERR $@,"\n";   
	return undef;
    }

    #print STDERR "starting $sf ",$sf->seqname,":",$sf->start,":",$sf->end,":",$sf->strand,"\n";
    
    
    # if this is something with subfeatures, then this is much more complex
    my @sub = $sf->sub_SeqFeature();

    #print STDERR "Got ",scalar(@sub),"sub features\n";

    if( $#sub >=  0 ) {
	# chain to constructor of the object. Not pretty this.
	$sf->flush_sub_SeqFeature();

	my $seen = 0;
	my $strand;
	foreach my $sub ( @sub ) {
	    #print STDOUT "Converting sub ",$sub->id,":",$sub->seqname,":",$sub->start,":",$sub->end,":",$sub->strand,"\n";
	    $sub = $self->_convert_seqfeature_to_vc_coords($sub);

	    if( !defined $sub ) {        
		next;
	    }
	    if( $seen == 0 ){
		$sf->start($sub->start);
		$sf->end($sub->end);
	    }

	    $seen =1;
	    $strand = $sub->strand;
	    $sf->add_sub_SeqFeature($sub,'EXPAND');
	}
	if( $seen == 1 ) {       
	    # we assumme that the mapping was unambiguous wrt to the strand

	    $sf->strand($strand);
	    
	    #print STDOUT "Giving back a new guy with start",$sf->start,":",$sf->end,":",$sf->strand," id ",$sf->id,"\n";
	    return $sf;
	} else {        
	    return undef;
	}
    }

    # might be clipped left/right
    #print ("Leftmost " . $mc->leftmost . " " . $mc->orientation . " " . $mc->start_in . " " . $mc->end_in  . " " . $sf->start . " " . $sf->end . "\n");


    if ($mc->leftmost){
	
	if ( $mc->orientation == 1) {
	    
	    # If end < startincontig contig for orientation 1 
	    if ($sf->start < $mc->start_in) {  
		return undef;              
	    }
	    
	} else {
	    # If start > startincontig for orientation <> 1
	    if ($sf->end > $mc->start_in) {  
		return undef;              
	    }
	}
    }  elsif ($mc->rightmost_end){
	
	if ( $mc->orientation == 1) {
	    
	    if ($sf->end >  $mc->rightmost_end) {  
		return undef;              
	    }
	    
	} else {
            # If start > startincontig for orientation <> 1
	    if ($sf->start <  $mc->rightmost_end) {  
		return undef;              
	    }
        }
    }
    
    
    # Could be clipped on ANY contig

    if( $mc->orientation == 1 ) {
	if ($sf->start < $mc->start_in) {  
	    return undef;              
	}
	if ($sf->end >  $mc->end_in) {  
	    return undef;              
	}
    } else {
	if ($sf->end > $mc->start_in) {  
	    return undef;              
	}
	if ($sf->start <  $mc->end_in) {  
	    return undef;              
	}
    }
	

    my ($rstart,$rend,$rstrand) = $self->_convert_start_end_strand_vc($cid,$sf->start,$sf->end,$sf->strand);
    
    $sf->start ($rstart);
    $sf->end   ($rend);
    $sf->strand($rstrand);
    
    if( $sf->can('attach_seq') ) {
	if (!$self->noseq) {
	    $sf->attach_seq($self->primary_seq);
	}
    }
    
    $sf->seqname($self->id);
    return $sf;
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

    my $mc;
    eval {
	$mc=$self->_vmap->get_MapContig($contig);
    };
    if($@) {
	$self->throw("Attempting to map a sequence feature with [$contig] on a virtual contig with no $contig");
    }

    if( $mc->orientation == 1 ) {
       
        # ok - forward with respect to vc. Only need to add offset
	my $offset = $mc->start - $mc->start_in;
	$rstart = $start + $offset;
	$rend   = $end + $offset;
	$rstrand = $strand;
    } else {
	my $offset = $mc->start+ $mc->start_in;
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
   print $fh "Contig Map Dump: \n";
   
   foreach my $mc ($self->_vmap->get_all_MapContigs) {
       print $fh "Contig ".$mc->contig->id." starts:",$mc->start," ends:".$mc->end." start in contig ",$mc->start_in," end in contig ".$mc->end_in." orientation ",$mc->orientation,"\n";
   }
}


=head2 get_all_RawContigs

 Title   : get_all_RawContigs
 Usage   : $obj->get_all_RawContigs()
 Function: 
 Example : 
 Returns : array of raw contigs sorted by starting position in vc
 Args    : 


=cut

sub get_all_RawContigs {
    my ($self) = @_;

    my @contigs = ();

    if( !ref $self || ! $self->isa('Bio::EnsEMBL::DB::VirtualContigI') ) {
        $self->throw("Must supply a VirtualContig to get_all_RawContigs: Bailing out...");
    }
    
    return $self->_vmap->get_all_RawContigs;
}

=head2 get_rawcontig_by_position

 Title   : get_rawcontig_by_position
 Usage   : $obj->get_rawcontig_by_position($position)
 Function: 
 Example : 
 Returns : returns a raw contig object or undef on error
 Args    : 


=cut

sub get_rawcontig_by_position {

    my ($self, $pos) = @_;
    
    if( !ref $self || ! $self->isa('Bio::EnsEMBL::DB::VirtualContigI') ) {
        $self->throw("Must supply a VirtualContig to get_all_RawContigs: Bailing out...");
    }
    
    if ($pos < 1 || $pos > $self->length){
        $self->throw("get_rawcontig_by_position error: Position must be > 0 and < vc->length");
    }

    foreach my $mc ($self->_vmap->get_all_MapContigs) {
        if ($pos > $mc->end) {
            next;
        } else {
            return $mc->contig;
        }
    }

    return (undef);
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
   my ($self,$cid) = @_;

   my $mc=$self->_vmap->get_MapContig($cid);
   return $mc->start;
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

   my $mc=$self->_vmap->get_MapContig($cid);
   return $mc->end;   
}

=head2 ori_in_vc

 Title   : ori_in_vc
 Usage   : $vc->ori_in_vc('rawcontigid')
 Function:
 Example :
 Returns : 
 Args    :


=cut
sub ori_in_vc {
    my ($self,$cid) = @_;
    
    my $mc=$self->_vmap->get_MapContig($cid);
    return $mc->orientation;   
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

=head2 _vmap

 Title   : _vmap
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _vmap{
    my ($obj,$value) = @_;
    
    if( defined $value) {
	$obj->{'_vmap'} = $value;
    }
    return $obj->{'_vmap'};    
}

=head2 noseq

 Title   : noseq
 Usage   : $vc->noseq
 Function: If set to 1, vc does not attach a primary seq to its features (exons) 
 Example : $vc->noseq(1)
 Returns : nothing
 Args    : 1/0


=cut

sub noseq{
    my ($obj,$value) = @_;
    
    if( defined $value) {
	$obj->{'_noseq'} = $value;
    }
    return $obj->{'_noseq'};    
}

1;






