
#
# BioPerl module for Bio::EnsEMBL::Virtual::Contig.pm
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Virtual::Contig - View a region of a genome as continuous DNA with features

=head1 SYNOPSIS

    # get a raw contig in the genome somehow.
    $rc = $db->get_Contig('AP000012.00001');

    # make a virtual contig from this position
    $vc = Bio::EnsEMBL::Virtual::Contig->new (
					      -focuscontig => $rc,
					      -focusposition => 100,
					      -orientation => 1,
					      -left => 1000000,
					      -right => 100000
					      );

    # this is now a ContigI compliant object, so the following calls work

    # This will be a long line!
    print "Sequence is ",$vc->primary_seq->seq,"\n";

    # if one part of a gene is on this then the gene is returned
    $vc->get_all_Genes();
   
    # this gets genes, but with start/end points in this virtual contig
    $vc->get_all_VirtualGenes();

    # can get out sequence features
    $vc->get_all_RepeatFeatures();
    $vc->get_all_SimilarityFeatures();
    $vc->get_all_PredictionFeatures();

    # from a virtual contig you can build either translations or 
    # sub virtual contigs

    # moves the VC 100,000 bp to the left
    $new_vc = $vc->extend(-100000,-100000);

    # builds a smaller vc, focused at 200,000 in this VC coordinates with
    # 1000 on each side
    $windowed_vc = $vc->(200000,1000,1000);


=head1 DESCRIPTION

This object provides a key abstraction for Ensembl: it presents a
series of underlying overlapping and ordered fragments which make up a
genome as a single piece of DNA with features. All the sequence,
sequence features and genes which are stored in the underlying
coordinates on the fragments (called RawContigs in Ensembl
terminology) can be "viewed" on this object, with coordinates all making sense
internally (though of course, comparing things between virtualcontigs is
simply crazy)

NB: To display a possible misunderstanding here, the fragments or
RawContigs are the local assemblies made by, for example, Phrap of
clones. In the future, as things move to larger assemblies, the
RawContigs are likely to the unit used in the computational
procedures.

This object is where the Ensembl software provides programmers with a
crucial abstraction which saves most people an very frustrating
programming task

=head1 AUTHOR - Ewan Birney

Email birney@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Virtual::Contig;
use vars qw(@ISA);
use strict;
# ContigI uses RootI.
use Bio::EnsEMBL::DB::ContigI;
use Bio::EnsEMBL::Virtual::Map;
use Bio::EnsEMBL::Virtual::PrimarySeq;

@ISA = qw(Bio::EnsEMBL::DB::ContigI);

sub _make_datastructures {
    my $self = shift;

    # data structures for caching and coordinating
    # top_SeqFeature call with EMBL dumping (ie, repressing BLAST hits in EMBL)
    $self->{'_sf_cache'} = {};
    $self->{'_skip_feature'} = {};
    my $vmap=Bio::EnsEMBL::Virtual::Map->new();
    $self->_vmap($vmap);
    return $self;
}

sub new {
    my ($class,@args) = @_;
    
    my $self = {};
    bless $self,$class;
    $self->_make_datastructures();

    my ($focuscontig,$focusposition,$ori,$leftsize,$rightsize,$clone) = 
	$self->_rearrange([qw( 
			       FOCUSCONTIG 
			       FOCUSPOSITION 
			       ORI 
			       LEFT 
			       RIGHT 
			       CLONE)],@args);

    #Create a new VirtualMap holder object for MapContigs
    my $vmap=Bio::EnsEMBL::Virtual::Map->new();
    $self->_vmap($vmap);

    # perhaps this should go into the Vmap constructor?
    $self->_vmap->right_overhang(0);
    $self->_vmap->left_overhang(0);
    my $VC_UNIQUE_NUMBER;

    if( defined $clone ) {
	$self->throw("Have not delt with clone VC rewrite yet");
    } else {
	# paranoia
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
	      $self->_vmap->build_map($focuscontig,$focusposition,$ori,$leftsize,$rightsize);
      $self->_vmap->dbobj($focuscontig->dbobj);
      $VC_UNIQUE_NUMBER = $focuscontig->id.".$focusposition.$ori.$leftsize.$rightsize";
    }
    my $length=$leftsize+$rightsize;
    print STDERR "CHECK!!!Got length $length\n";
    $self->length($length);
    $self->_unique_number($VC_UNIQUE_NUMBER);

    return $self;
}

		 
		    
=head1 New Virtual::Contig functions.

=head2 primary_seq

 Title   : primary_seq
 Usage   : $seq = $contig->primary_seq();
 Function: Gets a Bio::EnsEMBL::DB::VirtualPrimarySeq object out from the contig
 Returns : Bio::EnsEMBL::VirtualPrimarySeq object
 Args    : none

=cut

sub primary_seq {
    my ($self) = @_;

    my $vseq = Bio::EnsEMBL::Virtual::PrimarySeq->new( -vmap =>$self->_vmap, 
						       -id =>$self->id						        );
    return $vseq;
}

=head2 extend

 Title   : extend
 Usage   : $new_vc = $vc->extend(100,100);
 Function: Make a new vc by extending an existing one
 Example :
 Returns : Bio::EnsEMBL::Virtual::Contig
 Args    :


=cut

sub extend {
    my ($self, $left, $right) = @_;
    
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

=head2 extend_maximally

 Title   : extend_maximally
 Usage   : $new_vc = $vc->extend_maximally();
 Function: Extends an existing vc as far as possible in both directions
 Example :
 Returns : Bio::EnsEMBL::Virtual::Contig
 Args    :


=cut

sub extend_maximally {
    my ($self) = @_;
    
    # based on an original idea by Ewan Birney. ;)
    my $nvc = $self->extend(-10000000000,10000000000);
    return $nvc;
}

=head2 extend_maximally_left

 Title   : extend_maximally_left
 Usage   : $new_vc = $vc->extend_maximally_left();
 Function: Extends an existing vc as far as possible to the left
 Example :
 Returns : Bio::EnsEMBL::Virtual::Contig
 Args    :


=cut

sub extend_maximally_left {
    my ($self) = @_;
    # based on an original idea by Ewan Birney. ;)
    my $nvc = $self->extend(-10000000000,0);
    return $nvc;
}


=head2 extend_maximally_right

 Title   : extend_maximally_right
 Usage   : $new_vc = $vc->extend_maximally_right();
 Function: Extends an existing vc as far as possible to the right
 Example :
 Returns : Bio::EnsEMBL::Virtual::Contig
 Args    :


=cut

sub extend_maximally_right {
    my ($self) = @_;
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
    
    return Bio::EnsEMBL::Virtual::Contig->new(-focuscontig   => $map_contig,
					        -focusposition => $f_position,
					        -ori           => $ori,
					        -left          => $left,
					        -right         => $right
					        );
}

=head1 Functions implementing the Bio::SeqI interface inherieted by ContigI

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

=head2 length
    
 Title   : length
 Usage   : 
 Function: Provides the length of the contig
 Example :
 Returns : 
 Args    :


=cut

sub length {
    my $obj = shift;
    
    if( @_ ) {
	my $value = shift;
	$obj->{'length'} = $value;
    }
    return $obj->{'length'};
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

=head1 Functions to get sequence features and genes

=head2 get_all_VirtualGenes
    
 Title   : get_all_VirtualGenes
 Usage   : foreach my $virtualgene ( $contig->get_all_VirtualGenes ) 
 Function: Gets all the genes on this VC as VirtualGene objects
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

    # potentially we could be asked to convert something
    # that wasn't on this VC at all, eg, an exon from a distant contig
    eval {
	$mc=$self->_vmap->get_MapContig($cid);
    };
    if ($@) { 
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
	$mc=$self->_vmap->get_MapContig_by_id($contig);
    };
    if($@) {
	$self->throw("Attempting to map a sequence feature with [$contig] on a virtual contig with no $contig\n$@\n");
    }

    if( $mc->orientation == 1 ) {
       
        # ok - forward with respect to vc. Only need to add offset
	my $offset = $mc->start - $mc->rawcontig_start;
	$rstart = $start + $offset;
	$rend   = $end + $offset;
	$rstrand = $strand;
    } else {

	# flip strand
	$rstrand = $strand * -1;
	
	# yup. A number of different off-by-one errors possible here

	$rstart = $mc->end   - ($end   - $mc->rawcontig_start);
	$rend   = $mc->end   - ($start - $mc->rawcontig_start);

    }
    return ($rstart,$rend,$rstrand);
}


=head1 Helper functions for Virtual::Contig

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
   
   foreach my $mc ($self->_vmap->each_MapContig) {
       print $fh "Contig ".$mc->contig->id." starts:",$mc->start," ends:".$mc->end." start in contig ",$mc->rawcontig_start," end in contig ".$mc->rawcontig_end." orientation ",$mc->orientation,"\n";
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

=head1 Convertible Virtual::Contig stuff

=head2 convert_Gene_to_raw_contig

 Title   : convert_Gene_to_raw_contig
 Usage   : $newgene = $cvc->convert_Gene_to_raw_contig($gene)
 Function: Converts a gene built on this VirtualContig back to being
           a gene built on raw contig positions, read to be written back 
           

           Internally this builds a copy of the genes,transcripts and translations
           in RC coordinate space, making heavy use of VirtualMap raw_contig_position
           function. Exons could be split into sub exons across boundaries.

 Example :
 Returns : 
 Args    :


=cut

sub convert_Gene_to_raw_contig {
   my ($self,$gene) = @_;

   if( !ref $gene || ! $gene->isa('Bio::EnsEMBL::Gene') ) {
       $self->throw("Got to write a gene, not a [$gene]");
   }

   # sanity check 

   $self->_sanity_check($gene);

   # we need to map the exons back into RC coordinates.
   
   # in some cases the mapping will cause us to get a large number
   # of sticky exons.

   # the basically means we have to clone all the objects. :(

   my $clonedgene = Bio::EnsEMBL::Gene->new();
   my %translation;

   $clonedgene->id($gene->id);
   $clonedgene->version($gene->version);
   $clonedgene->created($gene->created);
   $clonedgene->modified($gene->modified);
   foreach my $dbl ( $gene->each_DBLink() ) {
       $clonedgene->add_DBLink($dbl);
   }

   foreach my $trans ( $gene->each_Transcript ) {
       my $clonedtrans = Bio::EnsEMBL::Transcript->new();
       $clonedtrans->id($trans->id);
       $clonedtrans->version($trans->version);
       $clonedtrans->created($trans->created);
       $clonedtrans->modified($trans->modified);
       foreach my $dbl ( $trans->each_DBLink() ) {
	   $clonedtrans->add_DBLink($dbl);
       }

       $clonedgene->add_Transcript($clonedtrans);

       foreach my $exon ( $trans->each_Exon ) {
	   my @clonedexons = $self->_reverse_map_Exon($exon);
	   foreach my $ce ( @clonedexons ) {
	       $clonedtrans->add_Exon($ce);
	   }
	   
	   # translations
	   if( exists $translation{$trans->translation->id} ) {
	       $clonedtrans->translation($translation{$trans->translation->id});
	   } else {
	       my $trl = $trans->translation(); 
	       my $clonedtrl = Bio::EnsEMBL::Translation->new();
	       $clonedtrl->id($trl->id);
	       $clonedtrl->start_exon_id($trl->start_exon_id);
	       $clonedtrl->end_exon_id($trl->end_exon_id);
	       $clonedtrl->version($trl->version);

	       my ($srawcontig,$start,$sstrand) = $self->_vmap->raw_contig_position($trl->start,1);
	       $clonedtrl->start($start);
	       my ($erawcontig,$end,$estrand) = $self->_vmap->raw_contig_position($trl->end,1);
	       $clonedtrl->end($end);
	       
	       $translation{$trl->id} = $clonedtrl;
	       $clonedtrans->translation($clonedtrl);
	   }
       }
   }

   return $clonedgene;

}

=head2 _reverse_map_Exon

 Title   : _reverse_map_Exon
 Usage   : (@exons) = $self->_reverse_map_Exon($exon)
 Function: Makes exons in RawContig coordinates from exon in VC coordinates.
           Multiple Exons might be returned when the Exons are made sticky
           due to exon crossing clone boundaries.
 Example :
 Returns : 
 Args    :


=cut

sub _reverse_map_Exon{
   my ($self,$exon) = @_;

   if( !ref $exon || !$exon->isa('Bio::EnsEMBL::Exon') ) {
       $self->throw("Must supply reverse map an exon not an [$exon]");
   }

   my ($scontig,$start,$sstrand) = $self->_vmap->raw_contig_position($exon->start,$exon->strand);
   my ($econtig,$end,$estrand)   = $self->_vmap->raw_contig_position($exon->end  ,$exon->strand);

   if( !ref $scontig || !ref $econtig || !$scontig->isa('Bio::EnsEMBL::DB::RawContigI') || !$econtig->isa('Bio::EnsEMBL::DB::RawContigI') ) {
       $self->throw("Exon on vc ".$exon->id." [".$exon->start.":".$exon->end."] is unmappable to rawcontig positions, probably being in a gap. Can't write");
   }

  
   if( $scontig->id eq $econtig->id ) {
       if( $sstrand != $estrand ) {
	   $self->throw("Bad internal error. Exon mapped to same contig but different strands!");
       }

       my $rmexon = Bio::EnsEMBL::Exon->new();
       $rmexon->id($exon->id);
       $rmexon->created($exon->created);
       $rmexon->modified($exon->modified);
       $rmexon->version($exon->version);
       $rmexon->phase($exon->phase);
       $rmexon->sticky_rank(1);
       foreach my $se ( $exon->each_Supporting_Feature ) {
	   my ($secontig,$sestart,$sestrand) = $self->_vmap->raw_contig_position($se->start,$se->strand);
	   my ($sncontig,$seend,$snstrand) = $self->_vmap->raw_contig_position($se->start,$se->strand);
	   if( !ref $secontig || !ref $sncontig || $secontig->id ne $sncontig->id ) {
	       $self->warn("supporting evidence spanning contigs. Cannot write");
	       next;
	   }
	   if( $sestart < $seend ) {
	       $se->start($sestart);
	       $se->end($seend);
	   } else {
	       $se->start($sestart);
	       $se->end($seend);
	   }
	   $se->strand($sestrand);
	   $se->seqname($secontig->id);
	   if( $se->can('attach_seq') ) {
	       $se->attach_seq($secontig->primary_seq);
	   }

	   $rmexon->add_Supporting_Feature($se);
       }

       # we could test on strand changes. This just assummes everything works
       # as it says on the tin ;)
       if( $start < $end ) {
	   $rmexon->start($start);
	   $rmexon->end($end);
       } else {
	   $rmexon->start($end);
	   $rmexon->end($start);
       }
       $rmexon->strand($sstrand);
       $rmexon->contig_id($scontig->id);
       $rmexon->seqname($scontig->id);
       $rmexon->attach_seq($scontig->primary_seq);
       return ($rmexon);
   } else {
       # we are in the world of sticky-ness....


       my @mapcontigs = $self->_vmap->each_MapContig();

       # walk to find scontig
       my $found = 0;
       my $mc;
       while ( $mc = shift @mapcontigs ) { 
 
	   if( $mc->contig->id eq $scontig->id ) {
	       print STDERR "Unshifting ",$mc->contig->id,"\n";
	       unshift(@mapcontigs,$mc);
	       $found = 1;
	       last;
	   }
       }
       if( $found == 0 ) {
	   $self->throw("Internal error - unable to find map contig with this id");
       }


       my $vcstart = $exon->start;
       #print STDERR "Looking from exon-wise",$exon->start,":",$exon->end,"\n";

       # ok. Move from start towards end, after we hit end.
       my @exported_exons;
       my $sticky = 1;

       foreach my $c ( @mapcontigs ) {	   
	   my $vcend;
	   #print STDERR "***Looking at",$c->contig->id," - $vcstart...\n";

	   if( $c->contig->id eq $econtig->id ) {
	       # go to end position
	       #print STDERR "Going for end...",$econtig->id,"\n";
	       $vcend = $exon->end();
	   } else {
	       #print STDERR "Going for end of contig\n";
	       $vcend = $c->end();
	   }

	   #print STDERR "....Going to call with $vcstart:$vcend\n";
	   #$self->_dump_map(\*STDERR);

	   my ($iscontig,$istart,$isstrand) = $self->_vmap->raw_contig_position($vcstart,$exon->strand);
	   my ($iecontig,$iend,$iestrand)   = $self->_vmap->raw_contig_position($vcend  ,$exon->strand);
  
	   if( $iscontig->id ne $iecontig->id || $isstrand != $iestrand) {
	       $self->throw("Bad internal error. Sticky Exon mapped to different contig/strand for a correct contig placement ".$iscontig->id.":".$iecontig->id);
	   }

	   my $rmexon = Bio::EnsEMBL::Exon->new();
	   $rmexon->id($exon->id);
	   $rmexon->created($exon->created);
	   $rmexon->modified($exon->modified);
	   $rmexon->version($exon->version);
	   $rmexon->phase($exon->phase);
	   $rmexon->sticky_rank($sticky++);
	   $rmexon->attach_seq($c->contig->primary_seq);

	   # could have flipped around, in which case, flip again.
	   # Again we assumme everything works as advertised
	   if( $istart > $iend ) {
	       $rmexon->start($iend);
	       $rmexon->end($istart);
	   } else {
	       $rmexon->start($istart);
	       $rmexon->end($iend);
	   }

	   $rmexon->strand($isstrand);
	   $rmexon->contig_id($c->contig->id);
	   $rmexon->seqname($c->contig->id);
	   push(@exported_exons,$rmexon);

	   if( $c->contig->id eq $econtig->id ) {
	       last;
	   }
	   $vcstart = $vcend+1;
       }

       return @exported_exons;
   }
       
   $self->throw("Internal error. Should not reach here!");
}

=head2 _sanity_check

 Title   : _sanity_check
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _sanity_check{
   my ($self,$gene) = @_;
   my $error =0;
   my $message;
   if( !defined $gene->id ) {
       $error = 1;
       $message .= "Gene has no id;";
   }
   if( !defined $gene->version ) {
       $error = 1;
       $message .= "Gene has to have a version;";
   }
   foreach my $trans ( $gene->each_Transcript ) {
       if( !defined $trans->id ) {
	   $error = 1;
	   $message .= "Transcript has no id;";
       }
       if( !defined $trans->translation || !ref $trans->translation) {
	   $error = 1;
	   $message .= "Transcript has no translation;";
       } else {
	   if( !defined $trans->translation->id ) {
	       $error = 1;
	       $message .= "Translation has no id";
	   } 
	   if( !defined $trans->translation->start ) {
	       $error = 1;
	       $message .= "Translation has no start";
	   } 
	   if( !defined $trans->translation->start_exon_id ) {
	       $error = 1;
	       $message .= "Translation has no start exon id";
	   } 
	   if( !defined $trans->translation->end ) {
	       $error = 1;
	       $message .= "Translation has no end";
	   } 
	   if( !defined $trans->translation->end_exon_id ) {
	       $error = 1;
	       $message .= "Translation has no end exon id";
	   } 
       }
       foreach my $exon ( $trans->each_Exon ) {
	   if( !defined $exon->id ) {
	       $error = 1;
	       $message .= "Exon has no id";
	   } 
	   if( !defined $exon->created ) {
	       $error = 1;
	       $message .= "Exon has no id";
	   } 
	   if( !defined $exon->modified ) {
	       $error = 1;
	       $message .= "Exon has no id";
	   } 
	   if( !defined $exon->contig_id  ) {
	       $error = 1;
	       $message .= "Exon has no contig id";
	   } else {
	       if( $exon->contig_id ne $self->id ) {
		   $error = 1;
		   $message .= "Exon [".$exon->id."] does not seem to be on this VirtualContig";
	       }
	   }
	   if( !defined $exon->start || !defined $exon->end || !defined $exon->strand || !defined $exon->phase ) {
	       $error = 1;
	       $message .= "Exon has error in start/end/strand/phase";
	   } 
       }
   }

   if( $error == 1 ) {
       $self->throw("Cannot write gene due to: $message");
   }
	   

}


=head1 Get/Set for Virtual::Contig orientated datastructures.

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

=head2 _unique_number

 Title   : _unique_number
 Usage   : $obj->_unique_number($newval)
 Function: 
 Example : 
 Returns : value of _unique_number
 Args    : newvalue (optional)


=cut

sub _unique_number{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'_unique_number'} = $value;
    }
    return $obj->{'_unique_number'};

}

=head2 _vmap

 Title   : _vmap
 Usage   : $obj->_vmap($newval)
 Function: 
 Example : 
 Returns : value of _vmap
 Args    : newvalue (optional)


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





