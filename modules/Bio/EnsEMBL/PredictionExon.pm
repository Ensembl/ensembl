#
# EnsEMBL module for Bio::EnsEMBL::PredictionExon
#
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=pod 

=head1 NAME Bio::EnsEMBL::PredictionExon - A class representing an Exon from an 
ab initio prediction method

=head1 SYNOPSIS

    $ex = new Bio::EnsEMBL::Exon(-START     => 100,
                                 -END       => 200,
                                 -STRAND    => 1,
                                 -SLICE     => $slice,
                                 -DBID      => $dbID,
				 -P_VALUE   => 23.5,
				 -SCORE     => 99
                                 );

   #seq returns a Bio::Seq
   my $seq = $exon->seq->seq();

   #peptide only makes sense within transcript context
   my $pep = $exon->peptide($transcript)->seq();

   #normal feature operations can be performed:
   $exon = $exon->transform('clone');
   $exon->move($new_start, $new_end, $new_strand);
   print $exon->slice->seq_region_name();

=head1 DESCRIPTION

This is a class which represents an prediction exon which is part of a 
predcition transcript. See Bio::EnsEMBL:PredictionTranscript

=head1 CONTACT

Post questions to the EnsEMBL developer list: <ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal
methods are usually preceded with a_

=cut

package Bio::EnsEMBL::PredictionExon;
use vars qw(@ISA);
use strict;


use Bio::EnsEMBL::Utils::Exception qw( warning throw deprecate );
use Bio::EnsEMBL::Utils::Argument qw( rearrange );


@ISA = qw(Bio::EnsEMBL::Exon);


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

  my ( $p_value, $score ) = 
    rearrange( [ "P_VALUE", "SCORE" ], @_ );

  $self->{'p_value'} = $p_value;
  $self->{'score'} = $score;

  return $self;
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
    throw( "End_phase setting not supported" );
  }
  return (3- $self->phase() + $self->length()) % 3;
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

=cut

sub transfer {
  my $self  = shift;
  
  my $new_exon = Bio::EnsEMBL::Feature::transform( $self, @_ );
  return undef unless $new_exon;

  #dont want to share the same sequence cache
  delete $new_exon->{'_seq_cache'};

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
  throw( "This object doesnt have supporting_feature" );
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
  throw( "This object doesnt have supporting_feature" );
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
  throw( "This object doesnt have supporting_feature" );
}




1;
