# EnsEMBL module for Transcript
# Copyright EMBL-EBI/Sanger center 2002
#
# Cared for by Arne Stabenau
#

#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

PredictionTranscript

=head1 SYNOPSIS


=head1 DESCRIPTION

Container for single transcript ab initio gene prediction ala GenScan.
Is directly storable/retrievable in EnsEMBL using PredictionTranscript Adaptor.
 

Creation:
   
     my $tran = new Bio::EnsEMBL::PredictionTranscript();
     $tran->add_Exon( $exon );

     my $tran = new Bio::EnsEMBL::PredictionTranscript(@exons);

     The order of the exons has to be right, as PT cant judge how to sort them.
     ( no sort as in Bio::EnsEMBL::Transcript )

     PredictionTranscript is geared towards the partial retrieve case from db.
     Exons can be missing in the middle. For storage though its necessary to 
     have them all and in contig coord system. 

Manipulation:

     # Returns an array of Exon objects, might contain undef instead of exon.
     my @exons = @{$tran->get_all_Exons};  

     # Returns the peptide translation as string 
     my $pep   = $tran->translate;

     # phase padded Exons cdna sequence. Phases usually match.
     my $cdna = $trans->get_cdna() 


=head1 CONTACT

contact EnsEMBL dev for information

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::PredictionTranscript;
use vars qw(@ISA);
use strict;


use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Transcript;
use Bio::Seq;

@ISA = qw(Bio::EnsEMBL::Root Bio::EnsEMBL::TranscriptI);

=head2 new

  Arg [1-]  : Bio::EnsEMBL::Exon $exon
              Exons which make up this transcript in the right order.
              Optional, can be added as well with add_Exon.
  Function  : Creates a new PredictionTranscript
  Returntype: Bio::EnsEMBL::PredictionTranscript
  Exceptions: none
  Caller    : Adaptor, Genscan Runnable

=cut

sub new {
  my($class,@optional_exons) = @_;

  if( ref $class ) { 
      $class = ref $class;
  }

  my $self = {};
  bless $self,$class;

  
  # set stuff in self from @args
  foreach my $exon (@optional_exons) {
    if( ! defined $self->{'exons'} ) {
      $self->{'exons'} = [];
    }
    
    $self->add_Exon($exon);
  }
  
  return $self;
}



=head2 stable_id

  Arg [1]    : (optional) $stable_id 
  Example    : my $pt_id = $prediction_transcript->stable_id;
  Description: Retrieves the stable id fro this prediction transcript.  
               Prediction transcripts do not maintain stable ids as real 
               transcripts do - the id is constructed from the name of the
               contig the prediction transcript was pulled off of, and the
               start and end of the transcript on the contig.
               i.e. the stable id is:  "$contig_name.$contig_start.$contig_end"
  Returntype : string
  Exceptions : none
  Caller     : general, PredictionTranscriptAdaptor

=cut

sub stable_id {
  my ($self, $value) = @_;

  if($value) {
    $self->{'_stable_id'} = $value;
  }

  return $self->{'_stable_id'};
}



=head2 coding_start

  Arg [1]  :  The new coding start of this prediction transcript in slice 
              coords.
  Function  : Getter/Setter for the coding start of this transcript.
              Implemented to satisfy requirements of TranscriptI interface
              and so that it can be drawn as a Transcript.
  Returntype: scalar int
  Exceptions: none
  Caller    : GlyphSet_transcript

=cut

sub coding_start {
  my ($self, $arg) = @_;

  if(defined $arg) {
    $self->{'coding_start'} = $arg;
  } elsif(!defined $self->{'coding_start'}) {
    #if the coding start is not defined, use the start of the transcript
    return $self->start();
  }

  return $self->{'coding_start'};
}


=head2 coding_end

  Arg [1]  :  (optional) The new coding end of this prediction transcript 
              in slice coords.
  Function  : Getter/Setter for the coding end of this transcript.
              Implemented to satisfy requirements of TranscriptI interface
              and so that it can be drawn as a Transcript.
  Returntype: scalar int
  Exceptions: none
  Caller    : GlyphSet_transcript

=cut

sub coding_end {
  my ($self, $arg) = @_;

  if(defined $arg) {
    $self->{'coding_end'} = $arg;
  } elsif(!defined $self->{'coding_end'}) {
    #if the coding end is not defined, use the end of the transcript
    return $self->end();
  }
  return $self->{'coding_end'};
}


=head2 start

  Arg [1]  :  The new start of this prediction transcript in slice 
              coords.
  Function  : Getter/Setter for the start of this transcript.
              Implemented to satisfy requirements of TranscriptI interface
              and so that it can be drawn as a Transcript.
  Returntype: scalar int
  Exceptions: none
  Caller    : GlyphSet_transcript

=cut

sub start {
  my ($self, $arg) = @_;

  if(defined $arg) {
    $self->{'start'} = $arg;
  } 

  return $self->{'start'};
}

=head2 end

  Arg [1]  :  The new end of this prediction transcript in slice coords.
  Function  : Getter/Setter for the end of this transcript.
              Implemented to satisfy requirements of TranscriptI interface
              and so that it can be drawn as a Transcript.
  Returntype: scalar int
  Exceptions: none
  Caller    : GlyphSet_transcript

=cut

sub end {
  my ($self, $arg) = @_;

  if(defined $arg) {
    $self->{'end'} = $arg;
  }

  return $self->{'end'};
}


 

## Attribute section ##

sub analysis {
  my ( $self, $value ) = @_;
  ( defined $value ) && 
    ( $self->{'analysis'} = $value );
  return $self->{'analysis'};
}


sub dbID {
   my $self = shift;
   
   if( @_ ) {
      my $value = shift;
      $self->{'dbID'} = $value;
    }
    return $self->{'dbID'};

}

sub adaptor {
   my $self = shift;
   
   if( @_ ) {
      my $value = shift;
      $self->{'adaptor'} = $value;
    }
    return $self->{'adaptor'};
}

## end attribute section ##

=head2 add_Exon

  Arg  1    : Bio::EnsEMBL::Exon $exon
  Arg [2]   : int $exon_position
              Use it when you know you dont have exons in the
              beginning. Useful for the Adaptor when retrieving 
              partial PT..
  Function  : Adds given Exon to this prediction transcript. 
              It can be at arbitrary position in the array. Not filled lower
              positions in the exon list are set undef then. Counting starts at 1.
  Returntype: none
  Exceptions: if argument is not Bio::EnsEMBL::Exon
  Caller    : Pipeline runnable Genscan

=cut


sub add_Exon {
  my ($self, $exon, $position) = @_;
  if( defined $position ) {
    $self->{'exons'}[$position-1] = $exon;
  } else {
    push( @{$self->{'exons'}}, $exon );
  }

  if(defined $exon && (!defined $self->{'start'} ||
		       $exon->start() < $self->{'start'})) {
    $self->start($exon->start());
  }
  if(defined $exon && (!defined $self->{'end'} ||
		       $exon->end() > $self->{'end'})) {
    $self->end($exon->end());
  }
}



=head2 get_all_Exons

  Args      : none
  Function  : Returns all Exons currently in the PredictionTranscript
              in the order 5' to 3'. If this is a partial PredictionTranscript,
              elements of the list will be undef.
  Returntype: listref Bio::EnsEMBL::Exon
  Exceptions: none
  Caller    : self->get_cdna(),Web for display.

=cut

sub get_all_Exons {
   my ($self) = @_;
   
   return $self->{'exons'};
}



=head2 get_all_translateable_Exons

  Arg [1]    : none
  Example    : $exons = $self->get_all_translateable_Exons
  Description: Retreives the same value of get_all_Exons for this prediction
               transcript.  In a prediction transcript there is no UTR and
               thus all exons are entirely translateable.
  Returntype : listref of Bio::EnsEMBL::Exon
  Exceptions : none
  Caller     : Embl_Dumper

=cut

sub get_all_translateable_Exons {
  my $self = shift;

  return $self->get_all_Exons;
}

############################################################


=head2 sort

 Function: Sorts the exons by start coordinate
           Sorts forward for forward strand and reverse for reverse strand
           It refills $self->{'exons'} with the sorted exons
 Returns : none
 Args    : none

=cut

sub sort {
  my $self = shift;

  # Fetch all the exons
  my @exons = @{$self->get_all_Exons()};

  # Empty the exon holder
  $self->flush_Exon();

  # Now sort the exons and put back in the feature table
  my $strand = $exons[0]->strand;

  if ($strand == 1) {
    @exons = sort { $a->start <=> $b->start } @exons;
  } 
  elsif ($strand == -1) {
    @exons = sort { $b->start <=> $a->start } @exons;
  }

  foreach my $e (@exons) {
    $self->add_Exon($e);
  }
}

############################################################

=head2 get_exon_count

  Args      : none
  Function  : How many exons are in this PTranscript. Some might 
              actually not be in this object, depending on how it was retrieved.
	      (non golden exons missing on Slice->get_predicitonTranscripts())     
  Returntype: int
  Exceptions: none
  Caller    : general

=cut


sub get_exon_count {
   my $self = shift;
   return scalar( @{$self->{'exons'}} );
}



=head2 set_exon_count

  Arg 1     : int $number_of_exons
              If the number of exons you put in with add_exon is not the 
              real number of exons in the Transcript, you can set it here.
              Might be necessary in db gets, where you dont get all.
  Function  : sets number of exons, so get_all_Exons returns more undef exons
              at the end. After this, you have to use position argument when you
	      want to insert Exons.
  Returntype: none
  Exceptions: If you set less then already in, you loose what you have :)
  Caller    : $self->adaptor()

=cut

sub set_exon_count {
  my $self = shift;
  my $number_of_exons = shift;

  if( ! defined $self->{'exons'} ) {
    $self->{'exons'} = [];
  }

  $#{$self->{'exons'}} = $number_of_exons-1;
}



=head2 length

  Args      : none
  Function  : length of DNA of all Exons in the PT together. No phase padding is done.
  Returntype: int
  Exceptions: differs from length in get_cdna() as that might be padded.
  Caller    : unknown

=cut


sub length {
    my( $self ) = @_;
    
    my $length = 0;
    foreach my $ex (@{$self->get_all_Exons}) {
      if( defined $ex ) { $length += $ex->length };
    }
    return $length;
}




=head2 flush_Exons

  Args      : none
  Function  : Removes all Exons from this PT. Resets count to zero.
  Returntype: none
  Exceptions: none
  Caller    : for completeness, use not likely.

=cut

sub flush_Exon{
   my ($self,@args) = @_;

   $self->{'exons'} = [];
   $self->{'start'} = undef;
   $self->{'end'} = undef;
}


=head2 translate

  Args      : none
  Function  : Give a peptide translation of all exons currently in
              the PT. Gives empty string when none is in.
  Returntype: a Bio::Seq as in transcript->translate()
  Exceptions: if the exons come in two or more groups, with an undef exon
              in the middle, only the first group is translated.
  Caller    : general

=cut


sub translate {
  my ($self) = @_;

  my $dna = $self->get_cdna();
  $dna    =~ s/TAG$|TGA$|TAA$//i;
  # the above line will remove the final stop codon from the mrna
  # sequence produced if it is present, this is so any peptide produced
  # won't have a terminal stop codon
  # if you want to have a terminal stop codon either comment this line out
  # or call translatable seq directly and produce a translation from it

  my $bioseq = new Bio::Seq( -seq => $dna, -moltype => 'dna' );
  
  return $bioseq->translate();
}
	 


=head2 get_cdna

  Args      : none
  Function  : Give a concat cdna of all exons currently in
              the PT. Gives empty string when none is in. Pads between not  
              phase matching exons. Builds internal coord translation table.
  Returntype: txt
  Exceptions: if the exons come in two or more groups, with an undef exon
              in the middle, only the first groups cdna is returned.
  Caller    : general, $self->translate()

=cut

sub get_cdna {
  my $self = shift;

  my $exons = $self->get_all_Exons();
  
  my $cdna = undef;
  my $lastphase = 0;

  my ( $cdna_start, $cdna_end );
  my ( $pep_start, $pep_end );
  my ( $new_cdna, $pep_count );

  $cdna_start = 1;
  $pep_start = 1;

  $self->{'_exon_align'} = [];

  for my $exon ( @$exons ) {
    my $exon_align = {};
    if( ! defined $exon ) {
      if( ! defined $cdna ) {
	next;
      } else {
	last;
      }
    } 

    push( @{$self->{'_exon_align'}}, $exon_align );

    my $phase = 0;

    if (defined($exon->phase)) {
      $phase = $exon->phase;
    }

    if( $phase != $lastphase ) {

      if( $lastphase == 1 ) {
	$cdna .= 'NN';
	$cdna_start += 2;
	$pep_start++;
      } elsif( $lastphase == 2 ) {
	$cdna .= 'N';
	$cdna_start += 1;
	$pep_start++;
      }

      #startpadding for this exon
      $cdna .= 'N' x $phase;
      $cdna_start += $phase;
    }
    
    $new_cdna = $exon->seq->seq();
#    print $new_cdna."\n";
    $cdna .= $new_cdna;
    $cdna_end = $cdna_start + CORE::length( $new_cdna ) - 1;

    # how many peptides are added by this exon??

    $pep_count = int( ( CORE::length( $new_cdna ) + $phase + 2 ) / 3 );

    $pep_end = $pep_start + $pep_count - 1; 
    $lastphase = $exon->end_phase();
      
    $exon_align->{ 'cdna_start' } = $cdna_start;
    $exon_align->{ 'cdna_end' } =  $cdna_end;
    $exon_align->{ 'pep_start' } = $pep_start;
    $exon_align->{ 'pep_end' } = $pep_end;
    $exon_align->{ 'exon' } = $exon;

    if( $lastphase == 0 ) { 
      $pep_start = $pep_end + 1;
    } else {
      $pep_start = $pep_end;
    }
    $cdna_start = $cdna_end+1;

  }

  if( ! defined $cdna ) { $cdna = '' };
  return $cdna
}



sub pep_coords {
    my $self = shift;
    $self->throw( "Not implemented yet" )
}


=head2 pep2genomic

  Arg  1    : int $start_amino
  Arg  2    : int $end_amino
  Function  : computes a list of genomic location covering the given range 
              of amino_acids. Exons do part of calculation. 
  Returntype: list [ start, end, strand, Contig_object, Exon, pep_start, pep_end ]
  Exceptions: none
  Caller    : Runnables mapping blast hits to genomic coords

=cut

sub pep2genomic {
   my $self = shift;
   my $start_amino = shift;
   my $end_amino = shift;
   my @result = ();

   my $start_cdna = $start_amino * 3 - 2;
   my $end_cdna = $end_amino * 3;

   my ( $ov_start, $ov_end );
   my ( $pep_start, $pep_end );

   print ::LOG "pep2genomic: ",$start_amino, " ", $end_amino,"\n";

   if( ! defined $self->{'_exon_align'} ) {
     $self->get_cdna();
   }

   for my $ex_align ( @{$self->{'_exon_align'}} ) {

     # calculate overlap of cdna region with each exon
     $ov_start = ( $start_cdna >= $ex_align->{cdna_start} ) ? $start_cdna : $ex_align->{cdna_start};
     $ov_end = ( $end_cdna <= $ex_align->{cdna_end} ) ? $end_cdna : $ex_align->{cdna_end};
     
     $pep_start = int(( $start_cdna+2)/3 );

     if( $ov_end >= $ov_start ) {
       my @ires = $ex_align->{exon}->cdna2genomic
	 ( $ov_start - $ex_align->{cdna_start} + 1,
	   $ov_end - $ex_align->{cdna_start} + 1 );
       for my $cdnaCoord ( @ires ) {
	 
	 push( @result, [ $cdnaCoord->[0], $cdnaCoord->[1], $cdnaCoord->[2],
			  $cdnaCoord->[3], $ex_align->{exon},
			  $cdnaCoord->[4] + $ex_align->{pep_start} - 1,
			  $cdnaCoord->[5] + $ex_align->{pep_start} - 1 ] );
       }
     }
   }
   for my $tmp ( @result ) {
     print ::LOG "pep2genomic result: ", join ( " ", @$tmp ),"\n";
   }

   return @result;
}


=head2 cdna2genomic

  Arg  1    : int $start_cdna
  Arg  2    : int $end_cdna
  Function  : computes a list of genomic location covering the given range 
              of nucleotides. Exons do part of calculation. 
  Returntype: list [ start, end, strand, Contig_object, Exon, cdna_start, cdna_end ]
  Exceptions: none
  Caller    : Runnables mapping blast hits to genomic coords

=cut

sub cdna2genomic {
   my $self = shift;
   my $start_cdna = shift;
   my $end_cdna = shift;
   my @result = ();


   my ( $ov_start, $ov_end );


   #print ::LOG "cdna2genomic: ",$start_cdna, " ", $end_cdna,"\n";

   if( ! defined $self->{'_exon_align'} ) {
     $self->get_cdna();
   }

   for my $ex_align ( @{$self->{'_exon_align'}} ) {

     # calculate overlap of cdna region with each exon
     $ov_start = ( $start_cdna >= $ex_align->{cdna_start} ) ? $start_cdna : $ex_align->{cdna_start};
     $ov_end = ( $end_cdna <= $ex_align->{cdna_end} ) ? $end_cdna : $ex_align->{cdna_end};
     
     if( $ov_end >= $ov_start ) {
       my @ires = $ex_align->{exon}->cdna2genomic
	 ( $ov_start - $ex_align->{cdna_start} + 1,
	   $ov_end - $ex_align->{cdna_start} + 1 );

       for my $cdnaCoord ( @ires ) {
	 
	 push( @result, [ $cdnaCoord->[0], $cdnaCoord->[1], $cdnaCoord->[2],
			  $cdnaCoord->[3], $ex_align->{exon},
			  $ov_start,
			  $ov_end ] );
       }
     }
   }
   for my $tmp ( @result ) {
     #print ::LOG "cdna2genomic result: ", join ( " ", @$tmp ),"\n";
   }

   return @result;
}


# debug helper
sub _dump {
  my $self = shift;
  my $res = "";
  if( ! defined $self->{'_exon_align'} ) {
    $self->get_cdna();
  }

  for my $ex ( @{$self->{'_exon_align'}} ) {
    $res .= "pep: ".$ex->{'pep_start'}." ".$ex->{pep_end}."\n";
    $res .= "exon: ".$ex->{'exon'}->start()." ".$ex->{'exon'}->end().
      " ".$ex->{'exon'}->strand()."\n";
  }
  return $res;
}

sub type {
  my ($self, $type) = @_;

  if(defined $type) {
    $self->{'_type'} = $type;
  }

  return $self->{'_type'};
}




1;
