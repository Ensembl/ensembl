#
# BioPerl module for Exon
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::StickyExon - This is a deprecated class that will be entirely
removed from the system. Do not use this class, use Bio::EnsEMBL::Exon instead

=head1 DESCRIPTION

StickyExons are deprecated. They should not be necessary and should not be
used. This module will be entirely removed from EnsEMBL in the near future.


=head1 CONTACT

The EnsEMBL developer mailing list for questions : <ensembl-dev@ebi.ac.uk>

=cut


package Bio::EnsEMBL::StickyExon;

use Bio::Seq;
use Bio::EnsEMBL::Utils::Exception qw(deprecate);

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::SeqFeature::Generic

use Bio::EnsEMBL::Exon;


@ISA = qw(Bio::EnsEMBL::Exon);


sub new {
  my($class,@args) = @_;
  
  my $self = Bio::EnsEMBL::Exon->new(@args);
  bless $self,$class;



  # Array to store exon tags
  $self->{_component_exons} = [];
  
  deprecate("StickyExon is a deprecated class.  It should not be needed." .
             "Use normal Exons instead.");

  return $self;
}


=head2 id

 Title   : id
 Usage   : overrides id to get/set locally and sets all component exons
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub id{
   my ($self,$value) = @_;

   if( defined $value ) {
       $self->{'_sticky_id'} = $value;
       foreach my $c ( @{$self->get_all_component_Exons} ) {
	   $c->id($value);
       }
   }

   return $self->{'_sticky_id'};

}



=head2 get_all_component_Exons

  Arg [1]    : 
  Example    : @component_exons = @{$exon->get_all_component_Exons}; 
  Description: Retrieves the component exons that comprise a sticky exon 
  Returntype : reference to a list of Bio::EnsEMBL::Exons
  Exceptions : none
  Caller     : general

=cut

sub get_all_component_Exons{
   my ($self,@args) = @_;

   return $self->{'_component_exons'};
}



=Head1 load_genomic_mapper

  Arg  1   : Bio::EnsEMBL::Mapper $mapper
             a mapper that will know hwo to go from cdna to genomic,
             after it is loaded here with the coordinates
  Arg  2   : int $id
             an id for the cdna, will probably be the address of the transcript
             that callewd this function. 

  Function : Loads the given mapper with cdna and genomic coordinates, so it can map 
             from one system to the other.

 Returntype: none
  Caller  : Bio::EnsEMBL::Transcript->convert_peptide_coordinate_to_contig


=cut


sub load_genomic_mapper {
  my ( $self, $mapper, $id, $start ) = @_;

  my $exons = $self->get_all_component_Exons;
  for my $exon ( @{$exons} ) {
    $mapper->add_map_coordinates( $id, $start, $start+$exon->length()-1,
				  $exon->strand(), $exon->contig,
				  $exon->start(), $exon->end() );
    $start += $exon->length;
  }
}

=head2 adjust_start_end

  Arg  1     : int $start_adjustment
  Arg  2     : int $end_adjustment
  Example    : none
  Description: returns a new Exon with this much shifted coordinates
  Returntype : Bio::EnsEMBL::Exon
  Exceptions : none
  Caller     : Transcript->get_all_translateable_Exons()

=cut


sub adjust_start_end {
  my ( $self, $start_adjust, $end_adjust ) = @_;

  if( $start_adjust == 0 && $end_adjust == 0 ) {
    return $self;
  }

  my $new_exon;
  
  my $start = 1 + $start_adjust;
  my $end = $self->length() + $end_adjust;

  my $mapper = Bio::EnsEMBL::Mapper->new( "cdna", "genomic" );
  my $current_start = 1;

  for my $exon ( @{$self->get_all_component_Exons()} ) {
    $mapper->add_map_coordinates( $self, $current_start, $current_start+$exon->length()-1,
                                  $exon->strand(), $exon->contig, $exon->start(), $exon->end() );
    $current_start += $exon->length();
  }

  my @mapped_coords = $mapper->map_coordinates( $self, $start, $end, 1, "cdna" );

  if( scalar @mapped_coords == 1 ) {
    # we can return a normal exon
    $new_exon = Bio::EnsEMBL::Exon->new();

    %$new_exon = %$self;
    $new_exon->start( $mapped_coords[0]->start() );
    $new_exon->end( $mapped_coords[0]->end() );
    $new_exon->strand( $mapped_coords[0]->strand() );
    $new_exon->contig( $mapped_coords[0]->id() );
    delete $new_exon->{'_component_exons'};    
  } else {
    # make a new sticky Exon
    $new_exon = Bio::EnsEMBL::StickyExon->new();
    %$new_exon = %$self;

    $new_exon->start( 1 );
    $new_exon->end( $end - $start + 1);
    $new_exon->strand( 1 );
    delete $new_exon->{'_component_exons'};    
    
    for my $coord ( @mapped_coords ) {
      my $cex = Bio::EnsEMBL::Exon->new();
      %$cex = %$self;
      $cex->start( $coord->start() );
      $cex->end( $coord->end() );
      $cex->strand( $coord->strand() );
      $cex->contig( $coord->id() ); 
      delete $cex->{'_component_exons'};
      $new_exon->add_component_Exon( $cex );
    }
  }

  return $new_exon;
}

=head2 add_component_Exon

 Title   : add_component_Exon
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_component_Exon{
   my ($self,$exon) = @_;

   if( !ref $exon || ! $exon->isa('Bio::EnsEMBL::Exon')) {
       $self->throw("$exon is not an exon");
   }

   push(@{$self->{'_component_exons'}},$exon);
}


sub contig {
  my ( $self ) = @_;

  return $self->{'_component_exons'}->[0]->contig();
}
  

=head2 length

 Title   : length
 Usage   : length
 Function: calculate number of  nucleotides constituting the Exon
 Example :
 Returns : a number
 Args    : none

=cut


sub length {
    my $self = shift;

    my $len =0; 

    foreach my $subexon ( @{$self->get_all_component_Exons} ) {
        $len += $subexon->length;
    }
    return $len;
}

=head2 _sort_by_sticky_rank

 Title   : _sort_by_sticky_rank
 Usage   : 
 Function: put the contained exons in the right order
 Example :
 Returns : 
 Args    : 

=cut


sub _sort_by_sticky_rank {
    my $self = shift;

    my @sorted;

    @sorted= sort {$a->sticky_rank <=> $b->sticky_rank } 
      @{$self->{'_component_exons'}};
    $self->{'_component_exons'} = \@sorted;
    return 1;
}




=head1 seq

  Arg [1]  : String $seq
             You can set the seq of a sticky Exon. 
             Omit, if you want to retrieve.
             contig method is probably the better way (sigh)
  Function : retrieve sequence of sticky exon from db or return stored one.
             If sequence was retrieved once, its cached.

  Returns  : Has to return Bio::Seq


=cut


sub seq {
  my $self = shift;
  my $seq = shift;

  if( defined $seq ) {
    if( $seq ) {
      $self->{'_seq'} = $seq;
    } else {
      $self->{'_seq'} = undef;
    }
  } else {
    my $seqString = "";
    for my $cExon ( @{$self->get_all_component_Exons} ) {
      $seqString .= $cExon->seq()->seq();
    }
    $self->{'_seq'} = $seqString;
  }

  return Bio::Seq->new( -seq => $self->{'_seq'} );
}



=head2 peptide

  Arg [1]    : Bio::EnsEMBL::Transcript $tr
  Example    : my $pep_str = $sticky_exon->peptide($transcript)->seq; 
  Description: StickyExon implementation of Bio::EnsEMBL::Exon::peptide.
               See Bio::EnsEMBL::Exon::peptide for details
  Returntype : Bio::Seq
  Exceptions : thrown if transcript argument is not provided
  Caller     : general

=cut

sub peptide {
  my $self = shift;
  my $tr   = shift;

  unless($tr && ref($tr) && $tr->isa('Bio::EnsEMBL::Transcript')) {
    $self->throw("transcript arg must be Bio::EnsEMBL:::Transcript not [$tr]");
  }
  
  my $pep_start = undef;
  my $pep_end   = undef;

  foreach my $exon (@{$self->get_all_component_Exons}) {
    #convert exons coordinates to peptide coordinates
    my @coords = 
      $tr->genomic2pep($exon->start, $exon->end, $exon->strand, $exon->contig);
    
    #filter out gaps
    @coords = grep {$_->isa('Bio::EnsEMBL::Mapper::Coordinate')} @coords;
 
    if(scalar(@coords) > 1) {
      $self->throw("Error. Exon maps to multiple locations in peptide." .
		   " Is this exon [$self] a member of this transcript [$tr]?");
      #if this is UTR then the peptide will be empty string
    } elsif(scalar(@coords) == 1) {
      my $c = $coords[0];
      #set the pep start to the minimum of all coords
      if(!defined $pep_start || $c->start < $pep_start) {
	$pep_start = $c->start;
      }

      #set the pep end to the maximum of all coords
      if(!defined $pep_end || $c->end > $pep_end) {
	$pep_end = $c->end;
      }
    }
  }

  #the peptide of this sticky is the region spanned by the component exons
  my $pep_str = '';
  if($pep_start && $pep_end) {
    $pep_str = $tr->translate->subseq($pep_start, $pep_end);
  }

  return Bio::Seq->new(-seq => $pep_str, 
		       -moltype => 'protein',
		       -alphabet => 'protein',
                       -id => $self->stable_id);
}



=head2 transform

  Arg  1    : Bio::EnsEMBL::Slice $slice
              make this slice coords.
  Function  : make slice coords from raw contig coords or vice versa
  Returntype: Bio::EnsEMBL::Exon (Bio::EnsEMBL::StickyExon)
  Exceptions: none
  Caller    : Gene::transform()

=cut


sub transform {
  my $self = shift;
  my $slice = shift;

  # print "Calling transform on sticky exon\n";
  
  if( defined $self->contig() and 
      $self->contig()->isa( "Bio::EnsEMBL::RawContig" ) )  {

    my $mapper = $slice->adaptor->db->get_AssemblyMapperAdaptor->fetch_by_type
      ( $slice->assembly_type() );

    my $dna_seq = "";
    my $mapped_start = 0;
    my $mapped_end = -1;
    my $composite_exon_strand = 0;
    my $composite_exon_phase;
    my $composite_exon_end_phase;

    my @supporting_features;
    
    # sort the component exons
    $self->_sort_by_sticky_rank(); 
    
    # and now retrieve them
    my $component_exons = $self->get_all_component_Exons;

    foreach my $c_exon ( @$component_exons ) {
     # print STDERR "component exon ".$self->stable_id." ".$c_exon->gffstring."\n";
      my @mapped = $mapper->map_coordinates_to_assembly
	(
	 $c_exon->contig()->dbID,
	 $c_exon->start(),
	 $c_exon->end(),
	 $c_exon->strand()
	);
      
      #    print STDERR "\nStickyExon.pm transform method:\n"; 
      #    print STDERR "[Mapped start  " . $mapped[0]->start . "\t";
      #    print STDERR "Mapped end    " . $mapped[0]->end . "\t";
      #    print STDERR "Mapped strand " . $mapped[0]->strand . "]\n";
      #    print STDERR "Sticky rank : " . $c_exon->sticky_rank() . "\n";
      
      # exons should always transform so in theory no error check
      # necessary
      if( ! @mapped ) {
	$self->throw( "Component Sticky Exon couldnt map" );
      }
    

      # should get a gap object returned if an exon lies outside of 
      # the current slice.  Simply return the exon as is - i.e. untransformed.
      # this untransformed exon will be distinguishable as it will still have
      # contig attached to it and not a slice.
      if( $mapped[0]->isa( "Bio::EnsEMBL::Mapper::Gap" )) {
	return $self;
      }

      # use empty Slice to map to chromosomal coords
      if( ! defined $slice->chr_name() ) {
	$slice->chr_name( $mapped[0]->id() );
	$slice->chr_start( 1 );
	$slice->strand( 1 );
      }

      # concatenate the raw sequence together

      # add the supporting features from the exons
      # each exon has the pieces of the supporting features that fall in the corresponding contig
      # they've been split before and at the moment they are not re-combined
      push @supporting_features, @{$c_exon->get_all_supporting_features()}; 
      
      #   print STDERR $c_exon->dbID . " " . $c_exon->seq . "\n";

      # now pull out the start and end points of the newly concatenated sequence
      # if we've got the first sticky exon, store the relevant info
      
      if ($c_exon->sticky_rank == 1 ) {
	# this assumes that the strand of the sticky exon is
	# set by the first one - is this correct?
	$composite_exon_strand = $mapped[0]->strand();
	
	if ( $composite_exon_strand == 1 ) {
	  # this means that in the forward strand the first component exon at the 5'end  in the slice
	  $mapped_start = $mapped[0]->start();
	  $composite_exon_phase = $c_exon->phase();
	}
	else { # for the reverse strand case it is something different
	  
	  # this means that in the reverse strand the first component exon at the 5'end  in the slice
	  $mapped_end = $mapped[0]->end();

	  # phase follows te 5' -> 3' sense, in contradistinction to start-end 
	  $composite_exon_phase = $c_exon->phase();
	}
      }
      
      # now do the end point
      # keep storing as you iterate over the component exons
      # since the exons are previously sorted based on their sticky rank
      # then the last set of stored values will be the last component exon
      else {
	if ( $mapped[0]->strand == 1 ) {	  
	  $mapped_end = $mapped[0]->end;
	  $composite_exon_end_phase  = $c_exon->end_phase;
	}
	else {
	  $mapped_start = $mapped[0]->start;
	  $composite_exon_end_phase = $c_exon->end_phase();
	}
      }
    }
    

    # now build the new composite exon
    my $newexon = Bio::EnsEMBL::Exon->new();

    if($slice->strand == 1) {
      $newexon->start ( $mapped_start - $slice->chr_start() + 1 );
      $newexon->end   ( $mapped_end   - $slice->chr_start() + 1);
      $newexon->strand( $composite_exon_strand );
    } 
    else {
      $newexon->start  ( $slice->chr_end() - $mapped_end   + 1);
      $newexon->end    ( $slice->chr_end() - $mapped_start + 1);
      $newexon->strand( $composite_exon_strand * -1);
    }
    $newexon->dbID($component_exons->[0]->dbID);
    $newexon->adaptor($component_exons->[0]->adaptor);

    $newexon->contig( $slice );
    
    #copy each of the supporting features and transform them
    my @feats;
    foreach my $sf (@supporting_features) {
      #my $f;
      #%$f = %$sf;
      #(mcvicker) this would be better as a copy
      push @feats, $sf->transform($slice);
    }
    $newexon->add_supporting_features(@feats);

    $newexon->phase( $composite_exon_phase );
    $newexon->end_phase( $composite_exon_end_phase );
    #print STDERR "transformed exon ".$newexon->stable_id." ".$newexon->gffstring."\n";
    return $newexon;
  } 
  else {
    $self->throw( "Unexpected StickyExon in Assembly coords ..." );
  }
}




=head2 get_all_supporting_features

  Arg [1]    : none
  Example    : @evidence = @{$sticky_exon->get_all_supporting_features};
  Description: Retreives any supporting features on this sticky exons 
               component exons. 
  Returntype : listreference of Bio::EnsEMBL::BaseAlignFeature objects 
  Exceptions : none
  Caller     : general

=cut

sub get_all_supporting_features {
  my $self = shift;

  my @out = ();

  foreach my $cexon (@{$self->get_all_component_Exons}) {
    push @out, @{$cexon->get_all_supporting_features};
  }

  return \@out;
}



=head2 add_supporting_features

  Arg [1]    : Bio::EnsEMBL::SeqFeatureI $feature
  Example    : $exon->add_supporting_features(@features);
  Description: Adds a list of supporting features to this exon. 
               Duplicate features are not added.  
               If supporting features are added manually in this
               way, prior to calling get_all_supporting_features then the
               get_all_supporting_features call will not retrieve supporting
               features from the database.
  Returntype : none
  Exceptions : throw if any of the features are not SeqFeatureIs
               throw if any of the features are not in the same coordinate
               system as the exon
  Caller     : general

=cut

sub add_supporting_features {
  my ($self,@features) = @_;

  # check whether this feature object has been added already
 FEATURE: foreach my $feature (@features) {
    foreach my $cexon (@{$self->get_all_component_Exons}) {
      
      if ((defined $cexon->contig() && defined $feature->contig())&&
	  ( $cexon->contig()->name() eq $feature->contig()->name())){
	$cexon->add_supporting_features($feature);
	next FEATURE;
      }
    }

    $self->warn("SupportingFeature could not be added, not on same contig " .
		"as component exons");
  }
}




1;
