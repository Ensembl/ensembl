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
. Is directly storable/retrievable in EnsEMBL using PredictionTranscript Adaptor. 

Creation:
   
     my $tran = new Bio::EnsEMBL::PredictionTranscript();
     $tran->add_Exon( $exon );

     my $tran = new Bio::EnsEMBL::PredictionTranscript(@exons);

     The order of the exons has to be right, as PT cant judge how to sort them.
     ( no sort as in Bio::EnsEMBL::Transcript )

     PredictionTranscript is geared towards the partial retrieve case from db. Exons can be missing in the middle. For storage though its necessary to have them all and in contig coord system. 

Manipulation:

     my @exons = $tran->get_all_Exons         # Returns an array of Exon objects, might contain undef instead of exon.
     my $pep   = $tran->translate()       # Returns the peptide translation as string
     my $cdna = $trans->get_cdna() # phase padded Exons cdna sequence.Padding might not often occur, phases usually match.


=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::PredictionTranscript;
use vars qw(@ISA);
use strict;


use Bio::Root::RootI;
use Bio::EnsEMBL::Exon;
use Bio::Seq;


@ISA = qw(Bio::Root::RootI);


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
  foreach my $a (@optional_exons) {
    if( ! defined $self->{'exons'} ) {
      $self->{'exons'} = [];
    }

    $self->add_Exon($a);
  }


  return $self;
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
}



=head2 get_all_Exons

  Args      : none
  Function  : Returns all Exons currently in the PredictionTranscript
              in the order 5' to 3'. If this is a partial PredictionTranscript,
              elements of the list will be undef.
  Returntype: list Bio::EnsEMBL::Exon
  Exceptions: none
  Caller    : self->get_cdna(),Web for display.

=cut

sub get_all_Exons {
   my ($self) = @_;
   
   return @{$self->{'exons'}};
}



=head2 get_exon_count

  Args      : none
  Function  : How many exons are in this PTranscript. Some might 
              actually not be in this object, depending on how it was retrieved.
	      (non golden exons missing on Slice->get_predicitonTranscripts())     Returntype: int
  Exceptions: none
  Caller    : general usage

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
    foreach my $ex ($self->get_all_Exons) {
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
}


=head2 translate

  Args      : none
  Function  : Give a peptide translation of all exons currently in
              the PT. Gives empty string when none is in.
  Returntype: txt
  Exceptions: if the exons come in two or more groups, with an undef exon
              in the middle, only the first group is translated.
  Caller    : general

=cut


sub translate {
  my ($self) = @_;

  my $dna = $self->get_cdna();
  my $bioseq = new Bio::Seq( -seq => $dna, -moltype => 'dna' );
  my $pep = $bioseq->translate();
  
  return $pep->seq();
}
	 


=head2 get_cdna

  Args      : none
  Function  : Give a concat cdna of all exons currently in
              the PT. Gives empty string when none is in. Pads between not  
              phase matching exons.
  Returntype: txt
  Exceptions: if the exons come in two or more groups, with an undef exon
              in the middle, only the first groups cdna is returned.
  Caller    : general, $self->translate()

=cut

sub get_cdna {
  my $self = shift;

  my @exons = $self->get_all_Exons();
  
  my $cdna = undef;
  my $lastphase = 0;

  for my $exon ( @exons ) {
    if( ! defined $exon ) {
      if( ! defined $cdna ) {
	next;
      } else {
	last;
      }
    } 
    
    if( $exon->phase() != $lastphase ) {

      if( $lastphase == 1 ) {
	$cdna .= 'NN';
      } elsif( $lastphase == 2 ) {
	$cdna .= 'N';
      }
      #startpadding for this exon
      $cdna .= 'N' x $exon->phase();
    }
    $cdna .= $exon->seq()->seq();
    $lastphase = $exon->end_phase();
  }

  if( ! defined $cdna ) { $cdna = '' };
  return $cdna
}


 
sub pep_coords {
    my $self = shift;
    $self->throw( "Not implemented yet" )
}


sub find_coord {
  my ($self,$coord,$type) = @_;
  $self->throw( "Not implemented yet" );
}



1;
