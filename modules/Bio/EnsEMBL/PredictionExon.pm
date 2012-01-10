=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::PredictionExon - A class representing an Exon from an ab
initio prediction method

=head1 SYNOPSIS

  $exon = new Bio::EnsEMBL::PredictionExon(
    -START   => 100,
    -END     => 200,
    -STRAND  => 1,
    -SLICE   => $slice,
    -DBID    => $dbID,
    -P_VALUE => 23.5,
    -SCORE   => 99
  );

  # seq() returns a Bio::Seq
  my $seq = $exon->seq->seq();

  # peptide() only makes sense within transcript context
  my $pep = $exon->peptide($transcript)->seq();

  # Normal feature operations can be performed:
  $exon = $exon->transform('clone');
  $exon->move( $new_start, $new_end, $new_strand );
  print $exon->slice->seq_region_name();

=head1 DESCRIPTION

This is a class which represents an prediction exon which is part of a
predcition transcript. See Bio::EnsEMBL:PredictionTranscript

=head1 METHODS

=cut

package Bio::EnsEMBL::PredictionExon;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Utils::Exception qw( warning throw deprecate );
use Bio::EnsEMBL::Utils::Argument qw( rearrange );


@ISA = qw(Bio::EnsEMBL::Exon);


=head2 new

  Args       : see SUPERCLASS Bio::EnsEMBL::Exon
  Example    : none
  Description: create an Exon object
  Returntype : Bio::EnsEMBL::PredictionExon 
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub new {
  my $class = shift;

  $class = ref $class || $class;

  my $self = $class->SUPER::new( @_ );

  my ( $p_value, $score ) = 
    rearrange( [ "P_VALUE", "SCORE" ], @_ );

  $self->{'p_value'} = $p_value;
  $self->{'score'} = $score;

  return $self;
}


=head2 score

  Arg [1]    : string $newval (optional) 
               The new value to set the score attribute to
  Example    : $score = $obj->score()
  Description: Getter/Setter for the score attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub score{
  my $self = shift;
  $self->{'score'} = shift if(@_);
  return $self->{'score'};
}



=head2 p_value

  Arg [1]    : string $newval (optional) 
               The new value to set the p_value attribute to
  Example    : $p_value = $obj->p_value()
  Description: Getter/Setter for the p_value attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub p_value{
  my $self = shift;
  $self->{'p_value'} = shift if(@_);
  return $self->{'p_value'};
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
  Status     : Stable

=cut



sub end_phase {
  my $self = shift;
  if( @_ ) {
    throw( "End_phase setting not supported" );
  }
  return ($self->phase() + $self->length()) % 3;
}


=head2 transform

  Arg  1     : String $coordinate_system_name
  Arg [2]    : String $coordinate_system_version
  Description: moves this exon to the given coordinate system. If this exon has
               attached supporting evidence, they move as well.
  Returntype : Bio::EnsEMBL::Exon
  Exceptions : wrong parameters
  Caller     : general
  Status     : Stable

=cut


sub transform {
  my $self = shift;

  # catch for old style transform calls
  if( !@_ || ( ref $_[0] && ($_[0]->isa( "Bio::EnsEMBL::Slice" ) or $_[0]->isa( "Bio::EnsEMBL::LRGSlice" )))) {
    throw( "transform needs coordinate systems details now," .
           "please use transfer" );
  }

  my $new_exon = Bio::EnsEMBL::Feature::transform( $self, @_ );
  return undef unless $new_exon;

  #dont want to share the same sequence cache
  delete $new_exon->{'_seq_cache'};

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
  Status     : Stable

=cut

sub transfer {
  my $self  = shift;

  my $new_exon = Bio::EnsEMBL::Feature::transfer( $self, @_ );
  return undef unless $new_exon;

  #dont want to share the same sequence cache
  delete $new_exon->{'_seq_cache'};

  return $new_exon;
}


=head2 add_supporting_features

  Description: For compatibility with Bio::EnsEMBL::Exon
               Does nothing.
  Returntype : none
  Status     : Stable

=cut

sub add_supporting_features { }


=head2 get_all_supporting_features

  Description: For compatibility with Bio::EnsEMBL::Exon
               Does nothing and returns empty list
  Returntype : empty list.
  Status     : Stable

=cut

sub get_all_supporting_features { return []; }


=head2 find_supporting_evidence

  Description: For compatibility with Bio::EnsEMBL::Exon
               Does nothing.
  Returntype : empty list.
  Status     : Stable

=cut

sub find_supporting_evidence { return []; }


1;
