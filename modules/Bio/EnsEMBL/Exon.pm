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

    $ex = new Bio::EnsEMBL::Exon(1208,1506,1) # start = 1208, end = 1506, forward strand
    
    Start and end coordinates are always stored with start < end.  If they are input
    in the reverse order they will be swapped over.  The value for the strand will
    be kept as its input value;

    Strand values:  + or  1 = forward strand
                    - or -1 = reverse strand
                    . or  0 = unknown strand

    $ex->attach_seq($dna);                    # $dna is a Bio::Seq
    $ex->phase(0);                            # Sets the phase of the exon
    $ex->end_phase();                         # Calculates the end_phase of the exon from the
                                              # Length of the dna and the starting phase

    Phase values  are 0,1,2

    $trans = $ex->translate();                # Translates the exon. Returns a Bio::Seq


    Frameshifts

    Frameshifts in the exon are stored as a multi-dimensional array of coordinates and lengths
    [start_position, length]

    my @fshifts =  $ex->get frameshifts();   # A multi-dim array of start coordinates and lengths

    Setting frameshifts
    
    $ex->add_frameshift(5,2);                # Adds a frameshift at position 5 and it
                                             # is an insertion of 2 bases
                                             # e.g CCCCAATTTT becomes cdna CCCCTTTT

    $ex->add_frameshift(5,-2);               # Adds a frameshift at position 5 and it
                                             # is a deletion of 2 bases
                                             # e.g CCCCAATTTT becomes cdna CCCCANNATTTT
                                             
    Retrieving cdna

    $dna = $ex->get_cdna();                  # Get a frameshift modified dna sequence
                                             # and return it in a Bio::Seq


=head1 DESCRIPTION

Exon object.  

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a_

=cut


# Let the code begin...


package Bio::EnsEMBL::Exon;
use vars qw(@ISA $AUTOLOAD);
use strict;

# Object preamble - inherits from Bio::SeqFeature::Generic

use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Slice;
use Bio::Seq; # exons have to have sequences...

@ISA = qw(Bio::EnsEMBL::SeqFeature);


sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);

  # Array to store supporting evidence for this exon

  # add in EnsEMBL tag as 1.

  #$self->primary_tag('exon');
  #$self->source_tag('EnsEMBL');

  # Parse the input paramters (start,end,strand)
  if ($#args == 2) {
      $self->_parse_args(@args);
  }
  # set exon rank to be 1 be default
  $self->sticky_rank(1);

  # set stuff in self from @args
  return $self; # success - we hope!
}

# Parse routine called from the constructor to
# set the basic variables start,end and strand

sub _parse_args {

  my ($self,$start,$end,$strand) = @_;

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

}

=pod 

=head1 Methods unique to exon


=pod 


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

  Arg [1]    : Bio::EnsEMBL::DBSQL::BaseAlignFeatureAdaptor $adaptor
  Example    : none
  Description: get/set for this objects Adaptor
  Returntype : Bio::EnsEMBL::DBSQL::BaseAlignFeatureAdaptor
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


=head2 contig_id

  Arg [1]    : string $contig_dbID
               You could set a contig ID in the Exon. Generally its a better
               idea to attach a complete contig to the exon.
  Example    : none
  Description: set a contig id (dbID or eventually others). The get part goes to the
               attached contig if there is nothing explicitly set.
  Returntype : string
  Exceptions : none
  Caller     : general

=cut


sub contig_id{
  my $self = shift;
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



=head2 contig

  Arg [1]    : Bio::EnsEMBL::RawContig/Slice $contig
  Example    : none
  Description: get/set for attribute contig. Is channeled to the bioperl
               comliant entire_seq/attach_seq calls.
  Returntype : Bio::EnsEMBL::RawContig/Slice
  Exceptions : none
  Caller     : general

=cut


sub contig {
  my $self = shift;
  if( @_ ) {
    my $value = shift;
    #print "setting exons contig to ".$value." \n";
    $self->attach_seq($value);
  }
  else {
    return $self->entire_seq();
  }
}


=head2 transform

  Arg  1    : Bio::EnsEMBL::Slice $slice
              make this slice coords
              if no slice, back to raw contig
  Function  : make slice coords from raw contig coords or vice versa
  Returntype: Bio::EnsEMBL::Exon (Bio::EnsEMBL::StickyExon)
  Exceptions: none
  Caller    : Gene::transform()

=cut

# could be used to transform from one slice to another ...
sub transform {
  my $self = shift;
  my $slice = shift;
  my $mapper;

  if( ! defined $slice ) {
    #Since slice arg is not defined -  we want raw contig coords
    if(( defined  $self->contig ) && 
       ( $self->contig->isa( "Bio::EnsEMBL::RawContig" )) ) {
      #we are already in rawcontig coords, nothing needs to be done
      return $self;
    } else {
      #transform to raw_contig coords from Slice coords
      return $self->_transform_to_rawcontig();
    }
  }

  #slice arg is defined - we want slice coords
  if( defined $self->contig ) {  
    if($self->contig->isa( "Bio::EnsEMBL::RawContig" ))  {
      #transform to slice coords from raw contig coords
      return $self->_transform_to_slice( $slice );
    } elsif($self->contig->isa( "Bio::EnsEMBL::Slice" )) {
      #transform to slice coords from other slice coords
      return $self->_transform_between_slices( $slice );
    } else {
      #Unknown contig type - throw an exception
      return $self->throw("Exon's 'contig' is of unknown type " 
		   . $self->contig() . " - cannot transform to Slice coords");
    }
  } else {
    #Can't convert to slice coords without a contig to work with
    return $self->throw("Exon's contig is not defined - cannot transform to " .
			"Slice coords");
  }
}



=head2 _transform_between_slices

  Arg  1     : Bio::EnsEMBL::Slice $new_slice
  Example    : none
  Description: Transforms the exons from one Slice to the given Slice, that needs to be
               on the same Chromosome
  Returntype : Bio::EnsEMBL::Exon
  Exceptions : Checks if Slice is attached and argument is Slice on same chromosome
  Caller     : transform

=cut



sub _transform_between_slices {
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
    #from slice is an empty slice, create a entire chromosome slice
    my $sa = $self->adaptor()->db()->get_SliceAdaptor();
    %$to_slice = %{$sa->fetch_by_chr_name($from_slice->chr_name())}

  }

  #sanity check - make sure we are transforming to the same chromosome
  if($to_slice->chr_name() ne $from_slice->chr_name()) {
    $self->throw("Cannot transform exon on chr " . $from_slice->chr_name() .
		 "to chr " . $to_slice->chr_name());
  }
  
  my $new_exon = new Bio::EnsEMBL::Exon();
  
  %$new_exon = %$self;

  my $shift = $from_slice->chr_start() - $to_slice->chr_start();

  #shift the start and end of the exon accordingly
  #negative coordinates are permitted, as are coords past slice boundary
  $new_exon->start($self->start() + $shift);
  $new_exon->end($self->end() + $shift);

  $new_exon->contig($to_slice);
  $new_exon->attach_seq($to_slice);

  return $new_exon;
}



=head2 _transform_to_slice

  Arg  1     : Bio::EnsEMBL::Slice $slice
  Example    : none
  Description: Transforms this Exon from RawContig coord to given Slice coord
  Returntype : Bio::EnsEMBL::Exon
  Exceptions : If the RawContig coords dont map
  Caller     : transform

=cut




sub _transform_to_slice {
  my $self = shift;
  my $slice = shift;

  
  my $mapper = $slice->adaptor->db->get_AssemblyMapperAdaptor->fetch_by_type
    ( $slice->assembly_type() );
  
  my @mapped = $mapper->map_coordinates_to_assembly
    (
     $self->contig()->dbID,
     $self->start(),
     $self->end(),
     $self->strand()
    );

  # exons should always transform so in theory no error check
  # necessary
  # actually we could have exons inside and outside the Slice because of db design
  # and the query that produces them
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
    my $sa = $self->adaptor()->db()->get_SliceAdaptor();
    %$slice = %{$sa->fetch_by_chr_name( $mapped[0]->id() )};
  } 

  my $newexon = Bio::EnsEMBL::Exon->new();
  %$newexon = %$self;
  
  $newexon->start( $mapped[0]->start() - $slice->chr_start() + 1);
  $newexon->end( $mapped[0]->end() - $slice->chr_start() + 1);
  $newexon->strand( $mapped[0]->strand() * $slice->strand() );
 

  $newexon->contig( $slice );
  $newexon->attach_seq( $slice );
  


  return $newexon;
}



=head2 _transform_to_rawcontig

  Args       : none
  Example    : none
  Description: Transform this Exon from Slice to RawContig coords
  Returntype : Bio::EnsEMBL::Exon
  Exceptions : Exon cant lie on Gap
  Caller     : transform

=cut



sub _transform_to_rawcontig {
  my $self = shift;

  my $mapper = $self->contig()->adaptor->db->get_AssemblyMapperAdaptor->fetch_by_type
    ( $self->contig()->assembly_type() );
  my $rcAdaptor = $self->adaptor()->db()->get_RawContigAdaptor();
  my $global_start = $self->contig->chr_start();

  
  $self->_transform_features_to_rawcontig();

  my @mapped = $mapper->map_coordinates_to_rawcontig
    (
     $self->contig()->chr_name(),
     $self->start()+$global_start-1,
     $self->end()+$global_start-1,
     $self->strand()*$self->contig()->strand()
    );

  if( ! @mapped ) {
    $self->throw( "Exon couldnt map" );
    return $self;
  }

  if( scalar( @mapped ) > 1 ) {
    # sticky exons
    # bjeacchh
   
    my $stickyExon = Bio::EnsEMBL::StickyExon->new();
    $stickyExon->phase( $self->phase() );
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
      $componentExon->contig_id( $rawContig->dbID );
      $componentExon->sticky_rank( $i + 1 );
      $componentExon->phase( $self->phase );
      $componentExon->end_phase($self->end_phase);
      $componentExon->dbID( $self->dbID() );
      $componentExon->adaptor( $self->adaptor() );
      $stickyExon->add_component_Exon( $componentExon );
      $sticky_length += ( $mapped[$i]->end() - $mapped[$i]->start() + 1 );
    }
    $stickyExon->end( $sticky_length );
    $stickyExon->strand( 1 );
    return $stickyExon;
    
  } else {
    # thats a simple exon
    if($mapped[0]->isa("Bio::EnsEMBL::Mapper::Gap")){
      
      $self->throw(" exon ". $self->start().
		   " ". $self->end(). " ".
		   $self->contig()->chr_name().
		   " lies on a gap cannot be mapped\n");
    }
    my $rawContig = $rcAdaptor->fetch_by_dbID( $mapped[0]->id() );
    $self->start( $mapped[0]->start() );
    $self->end( $mapped[0]->end() );
    $self->strand( $mapped[0]->strand() );
    # attaching seq ?
    $self->attach_seq( $rawContig );
    $self->contig( $rawContig );
    $self->contig_id($rawContig->dbID);
    return $self;
  }
}



# this should definatly happen inside the features ...
sub _transform_features_to_rawcontig {
  my ( $self ) = @_;

  if( ! defined $self->{'_supporting_evidence'} ) {
    return;
  }
  #print STDERR "TRANSFORMING SUPPORTING FEATURES\n";
  my @supporting_features = $self->each_Supporting_Feature;
  my @remapped_sf;
  foreach my $sf(@supporting_features) {
    #print "transforming ".$sf."\n";
    my @new_sf = $sf->transform();
    #print STDERR "have ".@new_sf." supporting features from transform\n";
    #foreach my $f(@new_sf){
    #  print "have ".$f." from transformation\n";
    #}
    push(@remapped_sf, @new_sf);
    
  }
  $self->{'_supporting_evidence'} = undef;
  foreach my $f(@remapped_sf){
   # print "adding ".$f." to supporting feature\n";
    $self->add_Supporting_Feature($f);
  }
  
}



=head2 pep_seq

  Title   : pep_seq
  Usage   : @pep = $feat->pep_seq
  Function: Returns the 3 frame peptide translations of the exon
  Returns : @Bio::Seq
  Args    : none

=cut

sub pep_seq {
    my ($self) = @_;
    my $pep = $self->_translate();
    return @$pep;
}

=head2 sticky_rank

 Title   : sticky_rank
 Usage   : $obj->sticky_rank($newval)
 Function: 
 Returns : value of sticky_rank
 Args    : newvalue (optional)


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

  Title   : end_phase
  Usage   : $end_phase = $feat->end_phase
  Function: Returns the end phase of the exon
  Returns : int
  Args    : none

=cut

sub end_phase {
  my ($self,$endphase) = @_;
  if ( defined($endphase) ){
    $self->{_end_phase} = $endphase;
  }
  if ( !defined( $self->{_end_phase} ) ){
    $self->warn("No end phase set in Exon. You must set it explicitly");
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

=cut

sub phase {
  my ($self,$value) = @_;
  
  if (defined($value)) {
    # Value must be 0,1,2, or -1 for non-coding
    if ($value =~ /^(-1|0|1|2)$/) {
#	print STDERR "Setting phase for " . $self->id . " to $value\n";
      $self->{'phase'} = $value;
    } else {
      $self->throw("Bad value ($value) for exon phase. Should only be -1,0,1,2\n");
    }
  }
  return $self->{'phase'};
}

=pod

=head2 frame

  Title   : frame
  Usage   : $frame = $feat->frame
  Function: Returns the frame of the exon in the parent DNA
  Returns : int
  Args    : none

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

=pod

=head2 type

  Title   : type
  Usage   : $type = $feat->type
  Function: Returns the type of the exon (Init,Intr,Term)
  Returns : String
  Args    : none

=cut

sub type {
  my ($self,$value) = @_;
  
  if (defined($value)) {
    $self->{'type'} = $value;
  }
  return $self->{'type'};
}


=pod

=head2 _translate

  Title   : _translate
  Usage   : @pep = $feat->_translate()
  Function: Returns the 3 frame translations of the exon
  Returns : @Bio::Seq
  Args    : none

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

=pod

=head2 translate

  Title   : translate
  Usage   : $pep = $feat->translate()
  Function: Returns the translation of the exon in the defined phase
  Returns : Bio::Seq
  Args    : none

=cut

sub translate {
    my($self) = @_;
    my $pep = $self->_translate() || throw ("Can't translate DNA\n");
    my $phase=$self->phase();
    if($phase){$phase=3-$phase;}
    return $pep->[$phase];
}

=head2 strand

 Title   : strand
 Usage   : $strand = $feat->strand()
 Function: Returns strand information, being 1,-1 or 0
 Returns : -1,1 or 0
 Args    : + or 1, - or -1, . or 0

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

 Title   : start_translation
 Usage   : $start_translation = $feat->start_translation
 Function: Returns coordinate taking into account phase
 Returns : number
 Args    : none

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


=head2 add_Supporting_Feature

 Title   : add_Supporting_Feature
 Usage   : $obj->add_Supporting_Feature($feature)
 Function: 
 Returns : Nothing
 Args    : Bio::EnsEMBL::SeqFeature


=cut

sub add_Supporting_Feature {
    my ($self,$feature) = @_;
    #print STDERR "args ".$feature."\n";
    $self->throw("Supporting evidence [$feature] not Bio::EnsEMBL::SeqFeatureI") unless 
	defined($feature) &&  $feature->isa("Bio::EnsEMBL::SeqFeatureI");
    #print STDERR "adding ".$feature." supporting feature\n";
    $self->{_supporting_evidence} = [] unless defined($self->{_supporting_evidence});

    # check whether this feature object has been added already
    my $found = 0;
    if ( $feature && $self->{_supporting_evidence} ){
      foreach my $added_feature ( @{ $self->{_supporting_evidence} } ){
	# compare objects
	if ( $feature == $added_feature ){
	  $found = 1;
	  
	  # no need to look further
	  last;
	}
      }
    }
    if ( $found == 0 ){
      #print STDERR "adding ".$feature." to evidence array\n";
      push(@{$self->{_supporting_evidence}},$feature);
    }
}


=head2 each_Supporting_Feature

 Title   : each_Supporting_Feature
 Usage   : my @f = $obj->each_Supporting_Feature
 Function: 
 Returns : @Bio::EnsEMBL::Feature
 Args    : none


=cut


sub each_Supporting_Feature {
    my ($self) = @_;
    if ( !defined ( $self->{_supporting_evidence} ) || 
	 (scalar @{$self->{_supporting_evidence}} == 0)) {
      $self->{_supporting_evidence} = [];
      if ( $self->adaptor ){
	$self->adaptor->fetch_evidence_by_Exon( $self );
      }
    }
    #print STDERR "returning ".@{$self->{_supporting_evidence}}." pieces of evidence\n";
    #foreach my $sf(@{$self->{_supporting_evidence}}){
    #  print STDERR "have ".$sf." evidence to return\n";
    #}
    return @{$self->{_supporting_evidence}};
}


=head2 find_supporting_evidence

 Title   : find_supporting_evidence
 Usage   : $obj->find_supporting_evidence(\@features)
 Function: Looks through all the similarity features and
           stores as supporting evidence any feature
           that overlaps with an exon.  I know it is
           a little crude but it\'s a start/
 Example : 
 Returns : Nothing
 Args    : Bio::EnsEMBL::Exon


=cut


sub find_supporting_evidence {
    my ($self,$features,$sorted) = @_;

    FEAT : foreach my $f (@$features) {
	# return if we have a sorted feature array
	if ($sorted == 1 && $f->start > $self->end) {
	    return;
	}
	if ($f->sub_SeqFeature) {
	  my @subf = $f->sub_SeqFeature;

	  $self->find_supporting_evidence(\@subf);
	} else {
	  if ($f->seqname eq $self->contig_id) {
	    if (!($f->end < $self->start || $f->start > $self->end)) {
	      $self->add_Supporting_Feature($f);
	    }
	  }
	}
      }
}


=head2 Ori methods

The ori methods are a hack around the virtual contig system to allow
us to cache genes in get_all_ExternalGenes via the ExternalWrapper
class, reusing the gene objects through more than one lift of the
virtual contig. Both Elia and Ewan understand this (vaguely) and it is
definitely a hack waiting to be removed somehow.


=head2 ori_start

 Title   : ori_start
 Usage   : $obj->ori_start($newval)
 Function: Getset for ori_start value
 Returns : value of ori_start
 Args    : newvalue (optional)


=cut

sub ori_start{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'ori_start'} = $value;
    }
    return $obj->{'ori_start'};

}

=head2 ori_end

 Title   : ori_end
 Usage   : $obj->ori_end($newval)
 Function: Getset for ori_end value
 Returns : value of ori_end
 Args    : newvalue (optional)


=cut

sub ori_end{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'ori_end'} = $value;
    }
    return $obj->{'ori_end'};

}

=head2 ori_strand

 Title   : ori_strand
 Usage   : $obj->ori_strand($newval)
 Function: Getset for ori_strand value
 Returns : value of ori_strand
 Args    : newvalue (optional)


=cut

sub ori_strand{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'ori_strand'} = $value;
    }
    return $obj->{'ori_strand'};

}


=head2 Stable id 

Stable id information is fetched on demand from stable tables

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
      my ($p,$f,$l) = caller;
      $self->warn("$f $l  created dates are loaded. Ignoring set value $value");
      return;
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
      my ($p,$f,$l) = caller;
      $self->warn("$f $l  modified dates are loaded. Ignoring set value $value");
      return;
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
      my ($p,$f,$l) = caller;
      $self->warn("$f $l  modified dates are loaded. Ignoring set value $value");
      return;
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
      return;
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
  
  my $phase_start_cdna = $start_cdna + $self->phase();
  my $phase_end_cdna = $end_cdna + $self->phase();

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


sub seq {
  my $self = shift;
  my $arg = shift;

  if( defined $arg ) {
    $self->{'_seq_cache'} = $arg;
  }

  if( defined $self->{'_seq_cache'} ) {
    return $self->{'_seq_cache'};
  }

  my $seq;
  #print STDERR " calling exon->seq\n";
  if ( ! defined $self->contig ) {
    $self->warn(" this exon doesn't have a contig you won't get a seq \n");
    return undef;
  }
  else {
    # call subseq on the contig which may be a RawContig or a Slice

#    print STDERR "[Exon.pm seq method: Start: " . $self->start . "\tEnd:   " . $self->end . "\t";
#    print STDERR "Strand: " . $self->strand . "]\nContig: " . $self->contig() . "\n\n";

      
    $seq = $self->contig()->subseq($self->start, $self->end);

    if($self->strand == -1){
      $seq =~ tr/ATGCatgc/TACGtacg/;
      $seq = reverse($seq);
    }
      
   }
  #print STDERR "have seq ".$seq."\n";
  my $bioseq = Bio::Seq->new(-seq=> $seq);
  $self->{'_seq_cache'} = $bioseq;

  return $bioseq;
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

=head2

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

=head2 id

 Title   : id
 Usage   : $obj->id($newval)
 Function: 
 Example : 
 Returns : value of id
 Args    : newvalue (optional)


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

1;
