
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
    
    return Bio::EnsEMBL::DB::VirtualContig->new(-focuscontig   => $map_contig,
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
    $self->throw('Ewan has not rewritten this function');

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
   $self->throw("Have not made this function...\n");
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





