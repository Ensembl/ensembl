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

use Bio::EnsEMBL::SeqFeature;
use Bio::Seq; # exons have to have sequences...
use Bio::EnsEMBL::StickyExon;

@ISA = qw(Bio::EnsEMBL::SeqFeature);



=head2 new

  Args       : see SUPERCLASS Bio::EnsEMBL::SeqFeature
  Example    : none
  Description: create an Exon object
  Returntype : Bio::EnsEMBL::Exon 
  Exceptions : none
  Caller     : general

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  
  # set exon rank to be 1 be default
  $self->sticky_rank(1);

  # set stuff in self from @args
  return $self; # success - we hope!
}



=head2 new_fast

  Arg [1]    : Bio::EnsEMBL::RawContig/Bio::EnsEMBL::Slice $contig
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
  my ($class,$contig,$start,$end,$strand) = @_;

  my $self = bless {}, $class;

  # Swap start and end if they're in the wrong order
  # We assume that the strand is correct and keep the input value.

  if ($start > $end) {
    my $tmp = $end;
    $end    = $start;
    $start  = $tmp;
  }
  
  $self->start ($start);
  $self->end   ($end);
  $self->strand($strand);
  $self->contig($contig);
  
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
      $self->{'adaptor'} = shift;
    }
    return $self->{'adaptor'};

}



=head2 _transform_between_Slices

  Arg [1]    : Bio::EnsEMBL::Slice $new_slice
  Example    : none
  Description: Transforms the exons from one Slice to the given Slice, 
               that needs to be on the same Chromosome. The method overwrites 
               the same method in Bio::EnsEMBL::SeqFeature
  Returntype : Bio::EnsEMBL::Exon
  Exceptions : Checks if Slice is attached and argument is Slice on same 
               chromosome.
  Caller     : transform

=cut

sub _transform_between_Slices {
  my ($self, $to_slice) = @_;

  my $from_slice = $self->contig();

  #sanity check - make sure we have something to transform from
  unless(defined $from_slice) {
    $self->throw("Exon's contig is not defined - cannot transform between "
		 . "slices\n");
  }
  #sanity check - make sure the from slice's chromosome is defined
  unless(defined $from_slice->chr_name()) {
    $self->throw("Exon's chromosome is not defined - cannot transform between "
		 . "slices\n");
  }

  unless(defined $to_slice->chr_name()) {
    #sanity check - we need an adaptor from a slice
    my $slice_adaptor = $to_slice->adaptor || $from_slice->adaptor;
    unless($slice_adaptor) {
      $self->throw("Exon cannot be transformed to empty slice without an " .
		   "an attached adaptor on the From slice or To slice");
    }
    #from slice is an empty slice, create a entire chromosome slice
    %$to_slice = %{$slice_adaptor->fetch_by_chr_name($from_slice->chr_name())};
  }

  #sanity check - make sure we are transforming to the same chromosome
  if($to_slice->chr_name() ne $from_slice->chr_name()) {
    $self->throw("Cannot transform exon on chr " . $from_slice->chr_name() .
		 "to chr " . $to_slice->chr_name());
  }
  
  #create a copy of the old exon
  my $new_exon = new Bio::EnsEMBL::Exon;
  %$new_exon = %$self;

  #unset the new exons supporting features
  delete $new_exon->{'_supporting_evidence'};

  #calculate the exons position in the assembly
  my ($exon_chr_start, $exon_chr_end, $exon_chr_strand);
  if($from_slice->strand == 1) {
    $exon_chr_start  = $self->start + $from_slice->chr_start - 1;
    $exon_chr_end    = $self->end   + $from_slice->chr_start - 1;
    $exon_chr_strand = $self->strand;
  } else {
    $exon_chr_start  = $from_slice->chr_end - $self->end   + 1;
    $exon_chr_end    = $from_slice->chr_end - $self->start + 1;
    $exon_chr_strand = $self->strand * -1;
  } 

  #now calculate the exons position on the new slice
  if($to_slice->strand == 1) {
    $new_exon->start($exon_chr_start - $to_slice->chr_start + 1);
    $new_exon->end  ($exon_chr_end   - $to_slice->chr_start + 1);
    $new_exon->strand($exon_chr_strand);
  } else {
    $new_exon->start($to_slice->chr_end - $exon_chr_end   + 1);
    $new_exon->end  ($to_slice->chr_end - $exon_chr_start + 1);
    $new_exon->strand($exon_chr_strand * -1);
  }

  $new_exon->contig($to_slice);

  #copy the attached supporting features and transform them
  my @feats;
  if( exists $self->{_supporting_evidence} ) {
    foreach my $sf (@{$self->get_all_supporting_features()}) {
      #my $f = $sf->new();
      #%$f = %$sf;
      ###(mcvicker) this would be better if the feature was copied
      push @feats, $sf->transform($to_slice);
    }
    $new_exon->add_supporting_features(@feats);
  }

  return $new_exon;
}



=head2 _transform_to_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
  Example    : none
  Description: Transforms this Exon from RawContig coord to given Slice coord. 
               The method overrides the same method in Bio::EnsEMBL::SeqFeature
  Returntype : Bio::EnsEMBL::Exon
  Exceptions : If the RawContig coords dont map
  Caller     : transform

=cut

sub _transform_to_Slice {
  my ($self,$slice) = @_;

  unless($self->contig) {
    $self->throw("Exon's contig must be defined to transform to Slice coords");
  }
  #print STDERR "transforming ".$self." from raw contig to slice coords\n";
  #print STDERR "exon ".$self->stable_id." ".$self->gffstring."\n";
  my $adaptor = $slice->adaptor || $self->contig->adaptor;

  unless($adaptor) {
    $self->throw("Cannot transform to exon slice unless either the " .
		 "exon->contig->adaptor or slice->adaptor is defined");
  }

  my $mapper = $adaptor->db->get_AssemblyMapperAdaptor->fetch_by_type
    ( $slice->assembly_type() );
  
  my @mapped = $mapper->map_coordinates_to_assembly
    (
     $self->contig()->dbID,
     $self->start(),
     $self->end(),
     $self->strand()
    );

  # exons should always transform so in theory no error check necessary
  # actually we could have exons inside and outside the Slice 
  # because of db design and the query that produces them
  if( @mapped != 1 ) {
    # Don't throw - return untransformed instead
    warn sprintf "Got %d features from map_coordinates_to_assembly() for exon %s",
        scalar(@mapped), $self->stable_id || $self->dbID;
    return $self;
  }
  # should get a gap object returned if an exon lies outside of 
  # the current slice.  Simply return the exon as is - i.e. untransformed.
  # this untransformed exon will be distinguishable as it will still have
  # contig attached to it and not a slice.
  if( $mapped[0]->isa( "Bio::EnsEMBL::Mapper::Gap" )) {
    return $self;
  }

  # the slice is an empty slice, create an enitre chromosome slice and
  # replace the empty slice with it
  if( ! defined $slice->chr_name() ) {
    my $slice_adaptor = $adaptor->db->get_SliceAdaptor;
    %$slice = %{$slice_adaptor->fetch_by_chr_name( $mapped[0]->id() )};
  } 

  my $newexon = new Bio::EnsEMBL::Exon();
  %$newexon = %$self;

  #unset supporting features of new exon
  delete $newexon->{'_supporting_evidence'};
  
  if ($slice->strand == 1) {
    $newexon->start( $mapped[0]->start() - $slice->chr_start() + 1);
    $newexon->end( $mapped[0]->end() - $slice->chr_start() + 1);
  } else {
    $newexon->start( $slice->chr_end() - $mapped[0]->end() + 1);
    $newexon->end( $slice->chr_end() - $mapped[0]->start() + 1);
  }

  $newexon->strand( $mapped[0]->strand() * $slice->strand() );
  $newexon->contig( $slice );
  #copy the attached supporting features and transform them
  my @feats;
  if( exists $self->{_supporting_evidence} ) {
    foreach my $sf (@{$self->get_all_supporting_features()}) {
      #my $f = $sf->new();
      #%$f = %$sf;
      #(mcvicker) this would be better if the feature was copied
      push @feats, $sf->transform($slice);
    }
    $newexon->add_supporting_features(@feats);
  }
  #print STDERR "transformed exon ".$newexon->stable_id." ".$newexon->gffstring."\n";
  return $newexon;
}



=head2 _transform_to_RawContig

  Args       : none
  Example    : none
  Description: Transform this Exon from Slice to RawContig coords
  Returntype : Bio::EnsEMBL::Exon
  Exceptions : Exon cant lie on Gap
  Caller     : transform

=cut

sub _transform_to_RawContig {
  my $self = shift;

  my $slice_adaptor = $self->contig->adaptor;

  unless($slice_adaptor) {
    $self->throw("Cannot transform exon to raw contig unless attached slice" .
		 " has adaptor defined. (i.e. exon->contig->adaptor)");
  }

  my $asma = $slice_adaptor->db->get_AssemblyMapperAdaptor();

  my $mapper = $asma->fetch_by_type( $self->contig()->assembly_type() );
  my $rcAdaptor       = $slice_adaptor->db->get_RawContigAdaptor();
  my $slice_chr_start = $self->contig->chr_start();
  my $slice_chr_end   = $self->contig->chr_end();

  my ($exon_chr_start,$exon_chr_end);

  if ($self->contig()->strand() == 1) {
    $exon_chr_start = $self->start() + $slice_chr_start - 1;
    $exon_chr_end   = $self->end()   + $slice_chr_start - 1;
  } 
  else {
    $exon_chr_end   = $slice_chr_end - $self->start() + 1;
    $exon_chr_start = $slice_chr_end - $self->end()   + 1;
  }

  my @mapped = $mapper->map_coordinates_to_rawcontig
    (
     $self->contig()->chr_name(),
     $exon_chr_start,
     $exon_chr_end,
     $self->strand()*$self->contig()->strand()
    );

  if( ! @mapped ) {
    $self->throw( "Exon couldnt map" );
    return $self;
  }

  #transform the supporting features to raw contig coords (hashed on contig)
  my %sf_hash;
  
  if( exists $self->{_supporting_evidence} ) {
    my $sfs = $self->get_all_supporting_features();
  SUPPORTING:foreach my $sf (@$sfs) {
      my @mapped_feats;
      eval{
	@mapped_feats = $sf->transform;
      };
      if($@){
	$self->warn("Supporting feature didn't mapped ignoring $@");
	next SUPPORTING;
      }
      foreach my $mapped_feat (@mapped_feats) {
	unless(exists $sf_hash{$mapped_feat->contig->name}) {
	  $sf_hash{$mapped_feat->contig->name} = [];
	}
	push @{$sf_hash{$mapped_feat->contig->name}}, $mapped_feat;
      }
    }
  }
      
  if( scalar( @mapped ) > 1 ) {
    # sticky exons
    # bjeacchh

    my $stickyExon = Bio::EnsEMBL::StickyExon->new();
    $stickyExon->phase( $self->phase() );
    $stickyExon->end_phase($self->end_phase());
    $stickyExon->adaptor( $self->adaptor() );
    $stickyExon->start( 1 );
    if( defined $self->dbID() ) { 
      $stickyExon->dbID( $self->dbID() ); 
    }

    my $sticky_length =0;
    # and then all the component exons ...
    for( my $i=0; $i <= $#mapped; $i++ ) {
      if($mapped[$i]->isa("Bio::EnsEMBL::Mapper::Gap")){
	$self->throw("exon '". $self->stable_id ."' lies on a gap cannot be mapped\n");
      }
      my $componentExon = Bio::EnsEMBL::Exon->new();
      $componentExon->start( $mapped[$i]->start() );
      $componentExon->end( $mapped[$i]->end() );
      $componentExon->strand( $mapped[$i]->strand() );
      my $rawContig = $rcAdaptor->fetch_by_dbID( $mapped[$i]->id() );
      $componentExon->contig( $rawContig );
      $componentExon->sticky_rank( $i + 1 );
      $componentExon->phase( $self->phase );
      $componentExon->end_phase($self->end_phase);
      $componentExon->dbID( $self->dbID() );
      $componentExon->adaptor( $self->adaptor() );

      #add the supporting features on this contig to the component exon
      if(exists $sf_hash{$rawContig->name}) {	
        $componentExon->add_supporting_features(@{$sf_hash{$rawContig->name}});
      }
      
      $stickyExon->add_component_Exon( $componentExon );
      $sticky_length += ( $mapped[$i]->end() - $mapped[$i]->start() + 1 );
    }
    $stickyExon->end( $sticky_length );
    $stickyExon->strand( 1 );
    if (defined($self->stable_id)) {
      $stickyExon->stable_id($self->stable_id); 
    }
    if (defined($self->version)) {
      $stickyExon->version($self->version);
    }
    if (defined($self->created)) {
      $stickyExon->created($self->created);
    }
    if (defined($self->modified)) {
      $stickyExon->modified($self->modified);
    }
    return $stickyExon;
    
  } else {
    # thats a simple exon
    if($mapped[0]->isa("Bio::EnsEMBL::Mapper::Gap")){
      
      $self->throw("exon '". $self->stable_id ."' lies on a gap cannot be mapped\n");
    }
    my $rawContig = $rcAdaptor->fetch_by_dbID( $mapped[0]->id() );
    my $new_exon = new Bio::EnsEMBL::Exon();
    
    #copy this exon
    %$new_exon = %$self;

    #unset supporting evidence of new exon
    delete $new_exon->{'_supporting_evidence'};

    $new_exon->start( $mapped[0]->start() );
    $new_exon->end( $mapped[0]->end() );
    $new_exon->strand( $mapped[0]->strand() );
    # attach raw contig
    $new_exon->contig( $rawContig );
    
    #replace old supporting feats with transformed supporting feats
    $new_exon->add_supporting_features(@{$sf_hash{$rawContig->name}});
    #print STDERR "transformed exon ".$new_exon->gffstring."\n";
    return $new_exon;
   
  }
}



=head2 sticky_rank

  Arg [1]    : (optional) int $value
  Example    : $sticky_rank = $exon->sticky_rank
  Description: Returns the position of a component exon in its assembled sticky
               exon.  Normal exons all have sticky_rank = 1 as do exons in
               slice coordinates. Exons in RawContig coordinates which span 
               multiple contigs are 'sticky' and are made of component exons
               whose order are denoted by this attribute.
  Returntype : int
  Exceptions : none
  Caller     : Bio::EnsEMBL::StickyExon

=cut

sub sticky_rank{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'sticky_rank'} = $value;
    }
    return $obj->{'sticky_rank'};

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
  my ($self,$endphase) = @_;
  if ( defined($endphase) ){
    $self->{_end_phase} = $endphase;
  }
  if ( !defined( $self->{_end_phase} ) ){
    $self->throw("No end phase set in Exon. You must set it explicitly. $!" .
	      "Caller: ".caller);
  }
  return $self->{_end_phase};
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

  $self->{_supporting_evidence} = [] 
    unless defined($self->{_supporting_evidence});
  
  # check whether this feature object has been added already
 FEATURE: foreach my $feature (@features) {
    #print STDERR "have ".$feature." to add to exon\n\n";
    unless($feature && $feature->isa("Bio::EnsEMBL::SeqFeatureI")) {
      $self->throw("Supporting feat [$feature] not a " . 
		   "Bio::EnsEMBL::SeqFeatureI");
    } 
    
    if ((defined $self->contig() && defined $feature->contig())&&
	    ( $self->contig()->name() ne $feature->contig()->name())){
      $self->throw("Supporting feat not in same coord system as exon\n" .
		   "exon is attached to [".$self->contig->name()."]\n" .
		   "feat is attached to [".$feature->contig->name()."]");
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




=head2 start

  Arg [1]    : int $start (optional)
  Example    : $start = $exon->start();
  Description: Getter/Setter for the start of this exon.  The superclass
               implmentation is overridden to flush the internal sequence
               cache if this value is altered
  Returntype : int
  Exceptions : none
  Caller     : general

=cut

sub start {
  my $self = shift;
  #if an arg was provided, flush the internal sequence cache
  delete $self->{'_seq_cache'} if(@_);
  return $self->SUPER::start(@_);
}


=head2 end

  Arg [1]    : int $end (optional)
  Example    : $end = $exon->end();
  Description: Getter/Setter for the end of this exon.  The superclass
               implmentation is overridden to flush the internal sequence
               cache if this value is altered
  Returntype : int
  Exceptions : none
  Caller     : general

=cut

sub end {
  my $self = shift;
  #if an arg was provided, flush the internal sequence cache
  delete $self->{'_seq_cache'} if(@_);
  return $self->SUPER::end(@_);
}


=head2 strand

  Arg [1]    : int $strand (optional)
  Example    : $start = $exon->strand();
  Description: Getter/Setter for the strand of this exon.  The superclass
               implmentation is overridden to flush the internal sequence
               cache if this value is altered
  Returntype : int
  Exceptions : none
  Caller     : general

=cut

sub strand {
  my $self = shift;
  #if an arg was provided, flush the internal sequence cache
  delete $self->{'_seq_cache'} if(@_);
  return $self->SUPER::strand(@_);
}

=head2 contig

  Arg [1]    : Bio::EnsEMBL::Slice or Bio::EnsEMBL::RawContig (optional)
  Example    : $slice = $exon->contig();
  Description: Getter/Setter for the contig of this exon.  The superclass
               implmentation is overridden to flush the internal sequence
               cache if this value is altered
  Returntype : int
  Exceptions : none
  Caller     : general

=cut

sub contig {
  my $self = shift;
  #if an arg was provided, flush the internal sequence cache
  delete $self->{'_seq_cache'} if(@_);
  return $self->SUPER::contig(@_);
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

    if(defined $value ) {
      $self->{'_created'} = $value;
    }


    if( exists $self->{'_created'} ) {
      return $self->{'_created'};
    }

    $self->_get_stable_entry_info();

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
    

    if( defined $value ) {
      $self->{'_modified'} = $value;
    }

    if( exists $self->{'_modified'} ) {
      return $self->{'_modified'};
    }

    $self->_get_stable_entry_info();

    return $self->{'_modified'};
}


=head2 version

 Title   : version
 Usage   : $obj->version()
 Function: 
 Returns : value of version
 Args    : 

=cut

sub version{

    my ($self,$value) = @_;
    

    if( defined $value ) {
      $self->{'_version'} = $value; 
    }

    if( exists $self->{'_version'} ) {
      return $self->{'_version'};
    }

    $self->_get_stable_entry_info();

    return $self->{'_version'};

}


=head2 stable_id

 Title   : stable_id
 Usage   : $obj->stable_id
 Function: 
 Returns : value of stable_id
 Args    : 


=cut

sub stable_id{

    my ($self,$value) = @_;
    

    if( defined $value ) {
      $self->{'_stable_id'} = $value;
    }

    if( exists $self->{'_stable_id'} ) {
      return $self->{'_stable_id'};
    }

    $self->_get_stable_entry_info();

    return $self->{'_stable_id'};

}

sub _get_stable_entry_info {
   my $self = shift;

   if( !defined $self->adaptor ) {
     return undef;
   }

   $self->adaptor->get_stable_entry_info($self);

}




=head1 load_genomic_mapper

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
    $tr->genomic2pep($self->start, $self->end, $self->strand, $self->contig);
  
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
    $pep_str = $pep->subseq($start, $end);
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

  if(!defined($self->{'_seq_cache'})) {
    my $seq;

    if ( ! defined $self->contig ) {
      $self->warn(" this exon doesn't have a contig you won't get a seq \n");
      return undef;
    }
    else {
      $seq = $self->contig()->subseq($self->start, $self->end);

      if($self->strand == -1){
        $seq =~ tr/ATGCatgc/TACGtacg/;
        $seq = reverse($seq);
      }
    }
    $self->{'_seq_cache'} = $seq;
  }

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



1;
