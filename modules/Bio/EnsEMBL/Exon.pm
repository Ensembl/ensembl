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

Bio::EnsEMBL::Exon - Confirmed Exon 

=head1 SYNOPSIS

    $ex = new Bio::EnsEMBL::Exon;

    $ex->start(10);
    $ex->end(100);

Examples of creating an exon

    # start = 1208, end = 1506, forward strand
    $ex = new Bio::EnsEMBL::Exon(1208,1506,1) 
    
    Start and end coordinates are always stored with start < end. If they are 
    input in the reverse order they will be swapped over.  The value for the 
    strand will be kept as its input value;

    Strand values:  + or  1 = forward strand
                    - or -1 = reverse strand
                    . or  0 = unknown strand

    $ex->contig($dna);     # $dna is a Bio::Seq
    $ex->phase(0);         # Sets the phase of the exon
    $ex->end_phase(1);      # sets the end_phase of the exon

    Phase values  are 0,1,2


=head1 DESCRIPTION

Exon object.  

=head1 CONTACT

Post questions to the EnsEMBL developer list: <ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal 
methods are usually preceded with a_

=cut


# Let the code begin...


package Bio::EnsEMBL::Exon;
use vars qw(@ISA $AUTOLOAD);
use strict;

# Object preamble - inherits from Bio::SeqFeature::Generic

use Bio::EnsEMBL::Feature;
use Bio::Seq; # exons have to have sequences...

use Bio::EnsEMBL::Utils::Exception qw( warning throw deprecate );
use Bio::EnsEMBL::Utils::Argument qw( rearrange );


@ISA = qw(Bio::EnsEMBL::Feature);


=head2 new

  Args       : see SUPERCLASS Bio::EnsEMBL::SeqFeature
  Example    : none
  Description: create an Exon object
  Returntype : Bio::EnsEMBL::Exon 
  Exceptions : none
  Caller     : general

=cut

sub new {
  my $class = shift;

  $class = ref $class || $class;

  my $self = $class->SUPER::new( @_ );
  
  my ( $phase, $end_phase, $stable_id, $version ) = 
    rearrange( [ "PHASE", "END_PHASE", "STABLE_ID", "VERSION" ], @_ );

  $self->{'phase'} = $phase;
  $self->{'end_phase'} = $end_phase;
  $self->{'stable_id'} = $stable_id;
  $self->{'version'} = $version;

  return $self;
}



=head2 new_fast

  Arg [1]    : Bio::EnsEMBL::Slice $slice
  Arg [2]    : int $start
  Arg [3]    : int $end
  Arg [4]    : int $strand (1 or -1)
  Example    : none
  Description: create an Exon object
  Returntype : Bio::EnsEMBL::Exon 
  Exceptions : none
  Caller     : general, creation in Bio::EnsEMBL::Lite::GeneAdaptor

=cut

sub new_fast {
  my ($class,$slice,$start,$end,$strand) = @_;

  my $self = bless {}, $class;

  # Swap start and end if they're in the wrong order
  # We assume that the strand is correct and keep the input value.

  if ($start > $end) {
    throw( "End smaller than start not allowed" );
  }
  
  $self->start ($start);
  $self->end   ($end);
  $self->strand($strand);
  $self->slice($slice);
  
  return $self;
}



=head2 dbID

  Arg [1]    : int $dbID
  Example    : none
  Description: get/set for the database internal id
  Returntype : int
  Exceptions : none
  Caller     : general, set from adaptor on store

=cut

sub dbID {
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'dbID'} = $value;
    }
    return $self->{'dbID'};

}



=head2 temporary_id

  Arg [1]    : string $temporary_id
  Example    : none
  Description: get/set for attribute temporary_id
               was invented from genebuild and shouldnt be necessary   
  Returntype : string
  Exceptions : none
  Caller     : general

=cut


sub temporary_id {
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'tempID'} = $value;
    }
    return $self->{'tempID'};

}


=head2 adaptor

  Arg [1]    : Bio::EnsEMBL::DBSQL::ExonAdaptor $adaptor
  Example    : none
  Description: get/set for this objects Adaptor
  Returntype : Bio::EnsEMBL::DBSQL::ExonAdaptor
  Exceptions : none
  Caller     : general, set from adaptor on store

=cut

sub adaptor {
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'adaptor'} = $value;
    }
    return $self->{'adaptor'};

}



=head2 end_phase

  Arg [1]    : (optional) int $end_phase
  Example    : $end_phase = $feat->end_phase;
  Description: Gets/Sets the end phase of the exon.
               end_phase = number of bases from the last incomplete codon of 
               this exon.
               Usually, end_phase = (phase + exon_length)%3
               but end_phase could be -1 if the exon is half-coding and its 3 
               prime end is UTR.
  Returntype : int
  Exceptions : warning if end_phase is called without an argument and the
               value is not set.
  Caller     : general

=cut



sub end_phase {
  my $self = shift;
  if( @_ ) { 
    $self->{'end_phase'} = shift;
  } else {
    if( ! defined ( $self->{'end_phase'} )) {
      warning( "No end phase set in Exon. You must set it explicitly. $!" );
    }
  }
  return $self->{'end_phase'};
}




=pod

=head2 phase

  my $phase = $exon->phase;
  $exon->phase(2);

Get or set the phase of the Exon, which tells the
translation machinery, which makes a peptide from
the DNA, where to start.

The Ensembl phase convention can be thought of as
"the number of bases of the first codon which are
on the previous exon".  It is therefore 0, 1 or 2
(or -1 if the exon is non-coding).  In ascii art,
with alternate codons represented by B<###> and
B<+++>:

       Previous Exon   Intron   This Exon
    ...-------------            -------------...

    5'                    Phase                3'
    ...#+++###+++###          0 +++###+++###+...
    ...+++###+++###+          1 ++###+++###++...
    ...++###+++###++          2 +###+++###+++...

Here is another explanation from Ewan:

Phase means the place where the intron lands
inside the codon - 0 between  codons, 1 between
the 1st and second base, 2 between the second and
3rd  base. Exons therefore have a start phase and
a end phase, but introns have just one phase.

=cut

sub phase {
  my ($self,$value) = @_;
  
  if (defined($value)) {
    # Value must be 0,1,2, or -1 for non-coding
    if ($value =~ /^(-1|0|1|2)$/) {
      #print STDERR "Setting phase to $value\n";
      $self->{'phase'} = $value;
    } else {
      $self->throw("Bad value ($value) for exon phase. Should only be" .
		   " -1,0,1,2\n");
    }
  }
  return $self->{'phase'};
}



=head2 frame

  Arg [1]    : none
  Example    : $frame = $exon->frame
  Description: Gets the frame of this exon
  Returntype : int
  Exceptions : thrown if an arg is passed
               thrown if frame cannot be calculated due to a bad phase value
  Caller     : general

=cut

sub frame {
  my ($self,$value) = @_;

  if( defined $value ) {
    $self->throw("Cannot set frame. Deduced from seq_start and phase");
  }

  # frame is mod 3 of the translation point

  if( $self->phase == -1 ) {
    return '.'; # gff convention for no frame info
  }
  if( $self->phase == 0 ) {
    return $self->start%3;
  }

  if( $self->phase == 1 ) {
    return ($self->start+2)%3;
  }

  if( $self->phase == 2 ) {
    return ($self->start+1)%3;
  }

  $self->throw("bad phase in exon ".$self->phase);

}



=head2 type

  Arg [1]    : (optional) $value
  Example    : Gets/Sets th etype of this exon
  Description: Returns the type of the exon (Init, Intr, Term)
  Returntype : string
  Exceptions : none
  Caller     : ?

=cut

sub type {
  my ($self,$value) = @_;
  
  if (defined($value)) {
    $self->{'type'} = $value;
  }
  return $self->{'type'};
}




=head2 transform

  Arg  1     : String $coordinate_system_name
  Arg [2]    : String $coordinate_system_version
  Description: moves this exon to the given coordinate system. If this exon has 
               attached supporting evidence, they move as well.
  Returntype : Bio::EnsEMBL::Exon
  Exceptions : wrong parameters
  Caller     : general

=cut


sub transform {
  my $self = shift;

  # catch for old style transform calls
  if( !@_ || ( ref $_[0] && $_[0]->isa( "Bio::EnsEMBL::Slice" ))) {
    throw( "transform needs coordinate systems details now, please use transfer" );
  }

  my $new_exon = $self->SUPER::transform( @_ );
  return undef unless $new_exon;

  if( exists $self->{'_supporting_evidence'} ) {
    my @new_features;
    for my $old_feature ( @{$self->{'_supporting_evidence'}} ) {
      my $new_feature = $old_feature->transform( @_ );
      push( @new_features, $new_feature );
    }
    $new_exon->{'_supporting_evidence'} = \@new_features;
  }
  return $new_exon;
}



=head2 transfer

  Arg [1]    : Bio::EnsEMBL::Slice $destination_slice
  Example    : none
  Description: Moves this Exon to given target slice coordinates. If Features
               are attached they are moved as well. Returns a new exon.
  Returntype : Bio::EnsEMBL::Gene
  Exceptions : none
  Caller     : general

=cut

sub transfer {
  my $self  = shift;
  
  my $new_exon = $self->SUPER::transfer( @_ );
  return undef unless $new_exon;

  if( exists $self->{'_supporting_evidence'} ) {
    my @new_features;
    for my $old_feature ( @{$self->{'_supporting_evidence'}} ) {
      my $new_feature = $old_feature->transfer( @_ );
      push( @new_features, $new_feature );
    }
    $new_exon->{'_supporting_evidence'} = \@new_features;
  }
  return $new_exon;
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
  #print STDERR "calling add supporting features\n\n";
  return unless @features;

  $self->{_supporting_evidence} ||= []; 
  
  # check whether this feature object has been added already
 FEATURE: foreach my $feature (@features) {
    #print STDERR "have ".$feature." to add to exon\n\n";
    unless($feature && $feature->isa("Bio::EnsEMBL::Feature")) {
      $self->throw("Supporting feat [$feature] not a " . 
		   "Bio::EnsEMBL::Feature");
    } 
    
    if ((defined $self->slice() && defined $feature->slice())&&
	    ( $self->slice()->name() ne $feature->slice()->name())){
      $self->throw("Supporting feat not in same coord system as exon\n" .
		   "exon is attached to [".$self->slice()->name()."]\n" .
		   "feat is attached to [".$feature->slice()->name()."]");
    }

    foreach my $added_feature ( @{ $self->{_supporting_evidence} } ){
      # compare objects
      if ( $feature == $added_feature ){
	#this feature has already been added
	next FEATURE;
      }
    }
    
    #no duplicate was found, add the feature
    push(@{$self->{_supporting_evidence}},$feature);
  }
}


=head2 get_all_supporting_features

  Arg [1]    : none
  Example    : @evidence = @{$exon->get_all_supporting_features()};
  Description: Retreives any supporting features added manually by 
               calls to add_supporting_features. If no features have been
               added manually and this exon is in a database (i.e. it h
  Returntype : listreference of Bio::EnsEMBL::BaseAlignFeature objects 
  Exceptions : none
  Caller     : general

=cut

sub get_all_supporting_features {
  my $self = shift;
  
  if( !exists  $self->{_supporting_evidence} )  {
    if($self->adaptor) {
      my $sfa = $self->adaptor->db->get_SupportingFeatureAdaptor();
      $self->{_supporting_evidence} = $sfa->fetch_all_by_Exon($self);
    } 
  }
  
  return $self->{_supporting_evidence} || [];
}


=head2 find_supporting_evidence

  Arg [1]    : Bio::EnsEMBL::SeqFeatureI $features
               The list of features to search for supporting (i.e. overlapping)
               evidence.
  Arg [2]    : (optional) boolean $sorted
               Used to speed up the calculation of overlapping features.  
               Should be set to true if the list of features is sorted in 
               ascending order on their start coordinates.
  Example    : $exon->find_supporting_evidence(\@features);
  Description: Looks through all the similarity features and
               stores as supporting features any feature
               that overlaps with an exon.  
  Returntype : none
  Exceptions : none
  Caller     : general

=cut

sub find_supporting_evidence {
  my ($self,$features,$sorted) = @_;

  foreach my $f (@$features) {
    # return if we have a sorted feature array
    if ($sorted == 1 && $f->start > $self->end) {
      return;
    }
    if ($f->sub_SeqFeature) {
      my @subf = $f->sub_SeqFeature;
      
      $self->find_supporting_evidence(\@subf);
    } 
    else {
      if ($f->entire_seq()->name eq $self->contig()->name) {
	if ($f->end >= $self->start && $f->start <= $self->end && $f->strand == $self->strand) {
	  $self->add_supporting_features($f);
	}
      }
    }
  }
}



=head2 created

 Title   : created
 Usage   : $obj->created()
 Function: 
 Returns : value of created
 Args    :


=cut

sub created{
    my ($self,$value) = @_;

    deprecated( "Created attribute not supported any more" );
    if(defined $value ) {
      $self->{'_created'} = $value;
    }


    return $self->{'_created'};

}

=head2 modified

 Title   : modified
 Usage   : $obj->modified()
 Function: 
 Returns : value of modified
 Args    : 


=cut

sub modified{
    my ($self,$value) = @_;
    

    deprecated( "Created attribute not supported any more" );
    if( defined $value ) {
      $self->{'_modified'} = $value;
    }


    return $self->{'_modified'};
}


=head2 stable_id

  Arg [1]    : string $stable_id
  Example    : none
  Description: get/set for attribute stable_id
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub stable_id {
   my $self = shift;
  $self->{'stable_id'} = shift if( @_ );
  return $self->{'stable_id'};
}


=head2 version

  Arg [1]    : string $version
  Example    : none
  Description: get/set for attribute version
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub version {
   my $self = shift;
  $self->{'version'} = shift if( @_ );
  return $self->{'version'};
}




=Head1 load_genomic_mapper

  Arg  1   : Bio::EnsEMBL::Mapper $mapper
             a mapper that will know how to go from cdna to genomic,
             after it is loaded here with the coordinates
  Arg  2   : int $id
             an id for the cdna, will probably be the address of the transcript
             that called this function. 

  Function : Loads the given mapper with cdna and genomic coordinates, so it can map 
             from one system to the other.

 Returntype: none
  Caller  : Bio::EnsEMBL::Transcript->convert_peptide_coordinate_to_contig


=cut


sub load_genomic_mapper {
  my ( $self, $mapper, $id, $start ) = @_;

  $mapper->add_map_coordinates( $id, $start, $start+$self->length()-1,
				$self->strand(), $self->contig,
				$self->start(), $self->end() );
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

  my $new_exon = Bio::EnsEMBL::Exon->new();
  %{$new_exon} = %{$self};

  #invalidate the sequence cache
  delete $new_exon->{'_seq_cache'};

  if( $self->strand() == 1 ) {
    $new_exon->start( $self->start() + $start_adjust );
    $new_exon->end( $self->end() + $end_adjust )
  } else {
    $new_exon->start( $self->start() - $end_adjust );
    $new_exon->end( $self->end() - $start_adjust )
  }

  return $new_exon;
}


=head2 peptide

  Arg [1]    : Bio::EnsEMBL::Transcript $tr
  Example    : my $pep_str = $exon->peptide($transcript)->seq; 
  Description: Retrieves the portion of the transcripts peptide
               encoded by this exon.  The transcript argument is necessary
               because outside of the context of a transcript it is not
               possible to correctly determine the translation.  Note that
               an entire amino acid will be present at the exon boundaries
               even if only a partial codon is present.  Therefore the 
               concatenation of all of the peptides of a transcripts exons 
               is not the same as a transcripts translation because the 
               summation may contain duplicated amino acids at splice sites.
               In the case that this exon is entirely UTR, a Bio::Seq object 
               with an empty sequence string is returned.
  Returntype : Bio::Seq
  Exceptions : thrown if transcript argument is not provided
  Caller     : general

=cut

sub peptide {
  my $self = shift;
  my $tr = shift;

  unless($tr && ref($tr) && $tr->isa('Bio::EnsEMBL::Transcript')) {
    $self->throw("transcript arg must be Bio::EnsEMBL:::Transcript not [$tr]");
  }

  #convert exons coordinates to peptide coordinates
  my @coords = 
    $tr->genomic2pep($self->start, $self->end, $self->strand, $self->slice);
  
  #filter out gaps
  @coords = grep {$_->isa('Bio::EnsEMBL::Mapper::Coordinate')} @coords;

  #if this is UTR then the peptide will be empty string
  my $pep_str = '';

  if(scalar(@coords) > 1) {
    $self->throw("Error. Exon maps to multiple locations in peptide." .
		 " Is this exon [$self] a member of this transcript [$tr]?");
  } elsif(scalar(@coords) == 1) {
    my $c = $coords[0];
    my $pep = $tr->translate;

    #bioperl doesn't give back residues for incomplete codons
    #make sure we don't subseq too far...
    my ($start, $end);
    $end = ($c->end > $pep->length) ? $pep->length : $c->end; 
    $start = ($c->start < $end) ? $c->start : $end;
    $pep_str = $tr->translate->subseq($start, $end);
  }
    
  return Bio::Seq->new(-seq => $pep_str, 
		       -moltype => 'protein',
		       -alphabet => 'protein',
                       -id => $self->stable_id);
}



=head2 seq

  Arg [1]    : none
  Example    : my $seq_str = $exon->seq->seq;
  Description: Retrieves the dna sequence of this Exon.  
               Returned in a Bio::Seq object.  Note that the sequence may
               include UTRs (or even be entirely UTR).
  Returntype : Bio::Seq
  Exceptions : warning if argument passed, warning if exon->contig not defined
  Caller     : general

=cut

sub seq {
  my $self = shift;
  my $arg = shift;

  if( defined $arg ) {
    $self->warn( "seq setting on Exon not supported currently" );
    $self->{'_seq_cache'} = $arg->seq();
  }

  if( defined $self->{'_seq_cache'} ) {
    return Bio::Seq->new(-seq=> $self->{'_seq_cache'});
  }

  my $seq;

  if ( ! defined $self->slice ) {
    $self->warn(" this exon doesn't have a slice you won't get a seq \n");
    return undef;
  }
  else {
      
    $seq = $self->slice()->subseq($self->start, $self->end);

    if($self->strand == -1){
      $seq =~ tr/ATGCatgc/TACGtacg/;
      $seq = reverse($seq);
    }
      
   }
  $self->{'_seq_cache'} = $seq;

  return Bio::Seq->new(-seq     => $self->{'_seq_cache'},
		       -id      => $self->stable_id,
		       -moltype => 'dna');
}


# Inherited methods
# but you do have all the SeqFeature documentation: reproduced here
# for convenience...

=pod

=head1 Methods inherited from SeqFeature

=head2 start

 Title   : start
 Usage   : $start = $feat->start
 Function: Returns the start coordinate of the feature
 Returns : integer
 Args    : none

=head2 end

 Title   : end
 Usage   : $end = $feat->end
 Function: Returns the end coordinate of the feature
 Returns : integer
 Args    : none

=head2 strand

 Title   : strand
 Usage   : $strand = $feat->strand()
           $feat->strand($strand)
 Function: get/set on strand information, being 1,-1 or 0
 Returns : -1,1 or 0
 Args    : none


=cut

=head2 length

 Title   : length
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=head2 sub_SeqFeature

 Title   : sub_SeqFeature
 Usage   : @feats = $feat->sub_SeqFeature();
 Function: Returns an array of sub Sequence Features
 Returns : An array
 Args    : none

=head2 primary_tag

 Title   : primary_tag
 Usage   : $tag = $feat->primary_tag()
 Function: Returns the primary tag for a feature,
           eg 'exon'
 Returns : a string 
 Args    : none

=head2 source_tag

 Title   : source_tag
 Usage   : $tag = $feat->source_tag()
 Function: Returns the source tag for a feature,
           eg, 'genscan' 
 Returns : a string 
 Args    : none

=head2 has_tag

 Title   : has_tag
 Usage   : $value = $self->has_tag('some_tag')
 Function: Returns the value of the tag (undef if 
           none)
 Returns : 
 Args    :

=head2 all_tags

 Title   : all_tags
 Usage   : @tags = $feat->all_tags()
 Function: gives all tags for this feature
 Returns : an array of strings
 Args    : none

=head2 gff_string

 Title   : gff_string
 Usage   : $str = $feat->gff_string
 Function: provides the feature information in GFF
           version 2 format.
 Returns : A string
 Args    : None



=head1 RangeI methods

These methods are inherited from RangeI and can be used
directly from a SeqFeatureI interface. Remember that a 
SeqFeature is-a RangeI, and so wherever you see RangeI you
can use a feature ($r in the below documentation).

=head2 overlaps

  Title   : overlaps
  Usage   : if($feat->overlaps($r)) { do stuff }
            if($feat->overlaps(200)) { do stuff }
  Function: tests if $feat overlaps $r
  Args    : a RangeI to test for overlap with, or a point
  Returns : true if the Range overlaps with the feature, false otherwise


=head2 contains

  Title   : contains
  Usage   : if($feat->contains($r) { do stuff }
  Function: tests whether $feat totally contains $r
  Args    : a RangeI to test for being contained
  Returns : true if the argument is totaly contained within this range


=head2 equals

  Title   : equals
  Usage   : if($feat->equals($r))
  Function: test whether $feat has the same start, end, strand as $r
  Args    : a RangeI to test for equality
  Returns : true if they are describing the same range


=head1 Geometrical methods

These methods do things to the geometry of ranges, and return
triplets (start, stop, strand) from which new ranges could be built.

=cut

=head2 intersection

  Title   : intersection
  Usage   : ($start, $stop, $strand) = $feat->intersection($r)
  Function: gives the range that is contained by both ranges
  Args    : a RangeI to compare this one to
  Returns : nothing if they don''t overlap, or 
            a new exon based on the range that they do overlap


=head2 union

  Title   : union
  Usage   : ($start, $stop, $strand) = $feat->union($r);
          : ($start, $stop, $strand) = Bio::RangeI->union(@ranges);
  Function: finds the minimal range that contains all of the ranges
  Args    : a range or list of ranges to find the union of
  Returns : the range containing all of the ranges

=cut

sub _get_stable_entry_info {
   my $self = shift;

   deprecated( "This function shouldnt be called any more" );
   if( !defined $self->adaptor ) {
     return undef;
   }

   $self->adaptor->get_stable_entry_info($self);

}



1;
