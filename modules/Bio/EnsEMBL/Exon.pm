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
    $ex->end_phase();      # Calculates the end_phase of the exon from the
                           # Length of the dna and the starting phase

    Phase values  are 0,1,2

    $trans = $ex->translate(); # Translates the exon. Returns a Bio::Seq


    Frameshifts

    Frameshifts in the exon are stored as a multi-dimensional array of 
    coordinates and lengths [start_position, length]

    # A multi-dim array of start coordinates and lengths
    my @fshifts =  $ex->get frameshifts();   

    Setting frameshifts
    
    $ex->add_frameshift(5,2);  # Adds a frameshift at position 5 and it
                               # is an insertion of 2 bases
                               # e.g CCCCAATTTT becomes cdna CCCCTTTT

    $ex->add_frameshift(5,-2); # Adds a frameshift at position 5 and it
                               # is a deletion of 2 bases
                               # e.g CCCCAATTTT becomes cdna CCCCANNATTTT
                                             
    Retrieving cdna

    $dna = $ex->get_cdna();    # Get a frameshift modified dna sequence
                               # and return it in a Bio::Seq


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
      my $value = shift;
      $self->{'adaptor'} = $value;
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
  $new_exon->{'_supporting_evidence'} = [];

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
  if( ! @mapped ) {
    $self->throw( "Exon couldnt map" );
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
  $newexon->{'_supporting_evidence'} = [];
  
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
    foreach my $sf (@$sfs) {
      foreach my $mapped_feat ($sf->transform()) {
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
	$self->throw(" exon lies on a gap cannot be mapped\n");
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
      
      $self->throw(" exon lies on a gap cannot be mapped\n");
    }
    my $rawContig = $rcAdaptor->fetch_by_dbID( $mapped[0]->id() );
    my $new_exon = new Bio::EnsEMBL::Exon();
    
    #copy this exon
    %$new_exon = %$self;

    #unset supporting evidence of new exon
    $new_exon->{'_supporting_evidence'} = undef;

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



=head2 pep_seq

  Arg [1]    : none
  Example    : @pep = $feat->pep_seq
  Description: Returns the 3 frame peptide translations of the exon
  Returntype : list of Bio::Seq objects
  Exceptions : none
  Caller     : ?

=cut

sub pep_seq {
    my ($self) = @_;
    my $pep = $self->_translate();
    return @$pep;
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



=head2 _translate

  Arg [1]    : none
  Example    : none
  Description: PRIVATE method Returns the 3 frame translations of the exon
  Returntype : listref of Bio::Seq objects
  Exceptions : thrown if exon length is less than 3
  Caller     : ?

=cut

sub _translate {
    my($self) = @_;
    my $pep;
    my $i;
  
    if( $self->length < 3 ) {
	$self->throw("Perfectly valid sub length 2 exon. Impossible to translate. Doh!".$self->length." ".$self->dbID);
    }

    # changed this to work with the new SeqFeature stuff. I am still not
    # 100% happy about this. EB.
  
    # Get the DNA sequence and create the sequence string
    $self->seq() || $self->throw("No DNA in object. Can't translate\n");

    my $dna = $self->seq();
  
    # Translate in all frames - have to chop
    # off bases from the beginning of the dna sequence 
    # for frames 1 and 2 as the translate() method
    # only translates in one frame. Pah!
  
    for ($i = 0; $i < 3; $i++) {
	my $tmp = new Bio::Seq(-seq => substr($dna->seq,$i));
	$pep->[$i] = $tmp->translate();
    }
    return $pep;
}



=head2 translate

  Arg [1]    : none
  Example    : $pep = $feat->translate;
  Description: Returns the translation of the exon in the defined phase
  Returntype : Bio::Seq
  Exceptions : thrown if the exon cannot be translated
  Caller     : ?

=cut

sub translate {
    my($self) = @_;

    my $pep = $self->_translate() || throw ("Can't translate DNA\n");

    my $phase= 0;

    if (defined($self->phase)) {
       $phase = $self->phase;
    }

    if ($phase){
       $phase = 3 - $phase;
    }

    return $pep->[$phase];
}



=head2 strand

  Arg [1]    : int $strand
  Example    : $strand = $exon->strand;
  Description: Gets/Sets strand of exon 
  Returntype : int
  Exceptions : none
  Caller     : general

=cut

sub strand {
  my ($self,$value) = @_;

  
   if( defined $value)  {
      if ($value eq "+") {$value =  1;}
      if ($value eq "-") {$value = -1;}
      if ($value eq ".") {$value =  0;}

      $self->{'strand'} = $value;
    }
    return $self->{'strand'};
}



=head2 start_translation

  Arg [1]    : none
  Example    : $start_translation = $exon->start_translation
  Description: Returns start translation coord taking into account phase
  Returntype : int
  Exceptions : thrown if arguments are passed
  Caller     : none

=cut

sub start_translation {
    my ($self,$value) = @_;
  
    if( defined $value){
	$self->throw("cannot set translation start!");
    }

    my $phase=$self->phase;
    $phase=3-$phase if $phase;
    if($self->strand==1){
	return ($self->start + $phase);
    }else{
	return ($self->end - $phase);
    }
}

=head2 end_translation

 Title   : end_translation
 Usage   : $end_translation = $feat->end_translation
 Function: Returns coordinate taking into account phase
 Returns : number
 Args    : none

=cut

sub end_translation {
    my ($self,$value) = @_;
  
    if( defined $value){
	$self->throw("cannot set translation end!");
    }

    my $phase=$self->end_phase;
    if($self->strand==1){
	return ($self->end - $phase);
    }else{
	return ($self->start + $phase);
    }
}

=head2 _rephase_exon_genscan

 Title   : _rephase_exon_genscan
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _rephase_exon_genscan{
   my ($self) = @_;

   my $dna = $self->seq();
   my $phase;
   my $pep = $self->_genscan_peptide();

   # trim top and bottem
   my $dnaseq = $dna->seq();

   LOOP : {

       $dna->setseq(substr($dnaseq,3,-3));
       
       my $tr1 = $dna->translate();
       my $tr1pep = $tr1->str();
       chop $tr1pep;

#        print STDERR "[$tr1pep] to\n[$pep]\n";
       if( $tr1pep !~ /\*/ && $pep =~ /$tr1pep/ ) {
	   $phase = 0;
	   last LOOP;
       }
       
       $dna->setseq(substr($dnaseq,4,-3));
       
       $tr1 = $dna->translate();
       $tr1pep = $tr1->str();
       chop $tr1pep;

#       print STDERR "[$tr1pep] to\n[$pep]\n";
       if( $tr1pep !~ /\*/ && $pep =~ /$tr1pep/ ) {
	   $phase = 1;
	   last LOOP;
       }
       
       $dna->setseq(substr($dnaseq,5,-3));
       
       $tr1 = $dna->translate();
       $tr1pep = $tr1->str();
       chop $tr1pep;

#       print STDERR "[$tr1pep] to\n[$pep]\n";
       if( $tr1pep !~ /\*/ && $pep =~ /$tr1pep/ ) {
	   $phase = 2;
	   last LOOP;
       }
       
   }

   

#   print STDERR "For exon ",$self->strand," ",$self->phase," ",$phase,"\n";
   $self->phase($phase);
   

}

=head2 _genscan_peptide

 Title   : _genscan_peptide
 Usage   : $obj->_genscan_peptide($newval)
 Function: 
 Returns : value of _genscan_peptide
 Args    : newvalue (optional)


=cut

sub _genscan_peptide{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_genscan_peptide'} = $value;
    }
    return $obj->{'_genscan_peptide'};

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
    else {
      $self->{_supporting_evidence} = [];
    }
  }
  
  return $self->{_supporting_evidence};
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


=head2 add_frameshift

 Title   : add_frameshift
 Usage   : $exon->add_frameshift($start,$length)
 Function: stores frameshift information in the current exon object
 Returns : Nothing
 Args    : start, length


=cut


sub add_frameshift {
    my ($self,$start,$length) = @_;


    # do some simple sanity checks

    my $tseq = $self->seq->seq;

    if ( $start > length( $tseq ) ) {
      print STDERR "Trying to add a frameshift outside the range of the sequence. Ignoring.\n";
      return;
    }

    if ( !defined $start ) {
      print STDERR "Trying to add a frameshift without any parameters. Ignoring.\n";
      return;
    }

    if ( !defined $length ) {
      print STDERR "Trying to add a frameshift without a specified length. Ignoring.\n";
      return;
    }

    if ( $start <= 0) {
      print STDERR "Trying to add a frameshift with a start point of zero or less. Ignoring.\n";
      return;
    }

    if ( $length == 0) {
      print STDERR "Trying to add a frameshift with a length of zero. Ignoring.\n";
      return;
    }

    # and finally if the frameshift data is valid....
    push @{$self->{'_frameshifts'}},[$start,$length];
}


=head2 get_frameshifts

 Title   : get_frameshifts
 Usage   : $exon->get_frameshifts()
 Function: retrieves an array of frameshifts if they exist
 Returns : Array of frameshifts [start, length] or UNDEF
 Args    : None

=cut

sub get_frameshifts {
    my ($self) = @_;

    # if already stored in the current exon object
    if ( defined $self->{'_frameshifts'} ) {
      return @{$self->{'_frameshifts'}};
      }
    else {  # else fetch frameshift data from database

      if ( !defined $self->adaptor ) {
	return;  # return undef if no ExonAdaptor
      }
      $self->adaptor->fetch_frameshifts( $self );
      if ( defined $self->{'_frameshifts'}){
	return @{$self->{'_frameshifts'}};
      }
      else {
	return;  # not necessarily any frameshifts stored
      }
    }
}


=head2 get_cdna

 Title   : get_cdna
 Usage   : $seq = $exon->get_cdna()
 Function: converts an exons dna based on any existing frameshifts
 Returns : Bio::Seq object
 Args    : None

=cut

sub get_cdna {
  my ($self) =@_;

  my $seq1 = $self->seq();
  my $seq = $seq1->seq();

  # check to see if frameshifts actually exist
  my @frameshifts = $self->get_frameshifts();

  # if the exon has no frameshifts, simply return the dna
  # sequence unmodified

  if ( scalar(@frameshifts) ==0  ) {
    my $temp_seq = Bio::Seq->new(
	   -SEQ         => $seq,
	   -DISPLAY_ID  => 'cdna_unmodified',
           -MOLTYPE     => 'dna'
           );

    return $temp_seq;
  }

  # if there are frameshifts...
  # sort frameshifts into order based on their start positions.
  # important that this is done because the while loop below relies on this!!
  # frameshifts are stored as [start, length]
  @frameshifts = sort { $a->[0] <=> $b->[0] } @frameshifts;

  my $fshift_seq = "";   # frameshift modified sequence

  # need to check that starting from 1 is in fact correct!!
  my $bp = 1;            # position along sequence
  my $curr_frame  = 0;    # current frameshift

  # run along the genomic sequence for the current exon
  while ( $bp <= length($seq) ) {

    # do something if we find a frameshift
    if ( $bp == $frameshifts[$curr_frame][0] ) {

      # if bps have been inserted then jump along the sequence
      # ignoring extraneous bps
      if ( $frameshifts[$curr_frame][1] > 0 ) {
	$bp += $frameshifts[$curr_frame][1];
      }

      # else if bps have been deleted, insert some Ns.
      # there shouldnt be a 0 case but cover it just in case.
      elsif ( $frameshifts[$curr_frame][1] <= 0 ) {
	$fshift_seq .= substr($seq, $bp-1, 1);
            # -1 to be consistent with the dna seqeunce	
      	$fshift_seq .= 'N' x -$frameshifts[$curr_frame][1];
	$bp++;
      }

      # to stop -warnings complaining
      if ( $curr_frame < $#frameshifts ) {
	$curr_frame++;      # point to the next frameshift
      }
   }
    else {  # store the current base pair
      $fshift_seq .= substr($seq, $bp-1, 1);
            # -1 to be consistent with the dna seqeunce
      $bp++;
    }
  }

  # a sanity check here to make sure the modified sequence is the
  # right length - but return modified cdna anyway...

  my $seq_changes = 0;

  for my $i ( 0 .. $#frameshifts ) {
    $seq_changes += $frameshifts[$i][1];
  }

  if ( length($seq) - $seq_changes != length($fshift_seq)) {
    print STDERR "Frameshift modified sequence isn't the correct length.\n";
  }

  # push the modified sequence into a Bio::Seq object
  my $temp_seq = Bio::Seq->new(
	 -SEQ         => $fshift_seq,
	 -DISPLAY_ID  => 'cdna_modified',
         -MOLTYPE     => 'dna'
         );

  return $temp_seq;
}



=head2 cdna2genomic

  Arg  1    : int $start_cdna
  Arg  2    : int $end_cdna
              relative to first nucleotide in this exon which is 1.
  Function  : calculates genomic coordinate for this range of cdna
              returns a list of [ $start, $end, $strand, $contig, start_pep, end_pep ]
  Returntype: list
  Exceptions: none
  Caller    : Transcript, PredictionTranscript

=cut

sub cdna2genomic {
  my $self = shift;
  my $start_cdna = shift;
  my $end_cdna = shift;
 
  my $phase = 0;

  if (defined($self->phase)) {
     $phase = $self->phase;
  } 
  my $phase_start_cdna = $start_cdna + $phase;
  my $phase_end_cdna = $end_cdna + $phase;

  my $pep_start = int(($phase_start_cdna+2)/3);
  my $pep_end = int (($phase_end_cdna+2)/3);

  if( $self->strand == 1 ) {
    return ([ $self->start + $start_cdna - 1,
	     $self->start + $end_cdna - 1,
	     $self->strand,
	     $self->contig,
	     $pep_start,
	     $pep_end ] );
  } else {
    return ( [ $self->end()- $end_cdna + 1,
	     $self->end() - $start_cdna + 1 ,
	     $self->strand(),
	     $self->contig,
	     $pep_start,  
	     $pep_end ] );
  }
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
    $tr->genomic2pep($self->start, $self->end, $self->strand, $self->contig);
  
  #filter out gaps
  my @coords = grep {$_->isa('Bio::EnsEMBL::Mapper::Coordinate')} @coords;

  #if this is UTR then the peptide will be empty string
  my $pep_str = '';

  if(scalar(@coords) > 1) {
    $self->throw("Error. Exon maps to multiple locations in peptide." .
		 " Is this exon [$self] a member of this transcript [$tr]?");
  } elsif(scalar(@coords) == 1) {
    my $c = $coords[0];
    $pep_str = $tr->translate->subseq($c->start, $c->end);
  }
    
  return Bio::Seq->new(-seq => $pep_str, 
		       -moltype => 'protein',
		       -alphabet => 'protein',
                       -id => $self->stable_id);
}



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

  return Bio::Seq->new(-seq=> $self->{'_seq_cache'});
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

=head2 length

 Title   : length
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

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



=head1 Deprecated Methods

=cut

###############################################################################

=head2 id

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use dbID or stable_id instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub id{
   my ($self) = shift;
   my $value  = shift;

   my ($p,$f,$l) = caller;
   $self->warn("$f:$l id deprecated. Please choose from stable_id or dbID");

  if( defined $value ) {
    $self->warn("$f:$l stable ids are loaded separately and dbIDs are generated on writing. Ignoring set value $value");
    return;
  }

   if( defined $self->stable_id ) {
     return $self->stable_id();
   } else {
     return $self->dbID;
   }
}



=head2 clone_id

  Args       : none
  Example    : none
  Description: deprecated, exons can have more than one clone.
               StickyExons dont support this call
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub clone_id{
   my $self = shift;

   $self->warn("Exon->clone_id method deprecated\n");
   return undef;

   if( @_ ) {
     my $value = shift;
     $self->{'clone_id'} = $value;
   }

   if( defined $self->{'clone_id'} ) {
     return $self->{'clone_id'};
   } elsif( defined $self->contig() ) {
     return $self->contig->cloneid();
   } else {
     return undef;
   }
}

=head2 contig_id

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use Bio::EnsEMBL::Exon::contig instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub contig_id{
  my $self = shift;
  $self->warn("Bio::EnsEMBL::Exon::contig_id is deprecated.  \n" .
	      "Use exon->contig->dbID instead $!");

#  if($contig_id) {
#    my $contig = 
#      $self->adaptor->db->get_RawContigAdaptor->fetch_by_dbID($contig_id);
#    $self->contig($contig);
#  }

#  return $self->contig()->dbID();
  
  if( @_ ) {
    my $value = shift;
    #print "setting contig_id = ".$value."\n";
    $self->{'contigid'} = $value;
  }
  if( defined $self->{'contigid'} ) {
    return $self->{'contigid'};
  } elsif( defined $self->contig() ) {
    return $self->contig->dbID();
  } else {
    return undef;
  }
}



=head2 each_Supporting_Feature

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use get_all_supporting_features instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub each_Supporting_Feature {
  my ($self, @args) = @_;
  
  $self->warn("Exon::each_Supporting_Feature has been renamed " .
	      "get_all_supporting_features" . caller);
  
  return $self->get_all_supporting_features(@args);
}


=head2 add_Supporting_Feature

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use add_supporting_features instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub add_Supporting_Feature {
  my ($self, @args) = @_;
  
  $self->warn("Exon::add_Supporting_Feature has been renamed " .
	      "add_supporting_features" . caller);

  return $self->add_supporting_features(@args);
}

=head2 ori_start

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED not needed, do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub ori_start{
   my $obj = shift;
   
   $obj->warn("Call to deprecated method ori_start " . caller);

   if( @_ ) {
      my $value = shift;
      $obj->{'ori_start'} = $value;
    }
    return $obj->{'ori_start'};

}


=head2 ori_end

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED not needed, do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub ori_end{
   my $obj = shift;

   $obj->warn("Call to deprecated method ori_end " . caller);

   if( @_ ) {
      my $value = shift;
      $obj->{'ori_end'} = $value;
    }
    return $obj->{'ori_end'};

}


=head2 ori_strand

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED not needed, do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub ori_strand{
   my $obj = shift;

   $obj->warn("Call to deprecated method ori_strand " . caller);

   if( @_ ) {
      my $value = shift;
      $obj->{'ori_strand'} = $value;
    }
    return $obj->{'ori_strand'};
}


1;
