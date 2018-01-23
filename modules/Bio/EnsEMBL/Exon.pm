=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Exon - A class representing an Exon

=head1 SYNOPSIS

    $exon = new Bio::EnsEMBL::Exon(
      -START     => 100,
      -END       => 200,
      -STRAND    => 1,
      -SLICE     => $slice,
      -DBID      => $dbID,
      -ANALYSIS  => $analysis,
      -STABLE_ID => 'ENSE000000123',
      -VERSION   => 2
    );

  # seq() returns a Bio::Seq
  my $seq = $exon->seq->seq();

  # Peptide only makes sense within transcript context
  my $pep = $exon->peptide($transcript)->seq();

  # Normal feature operations can be performed:
  $exon = $exon->transform('clone');
  $exon->move( $new_start, $new_end, $new_strand );
  print $exon->slice->seq_region_name();

=head1 DESCRIPTION

This is a class which represents an exon which is part of a transcript.
See Bio::EnsEMBL:Transcript

=head1 METHODS

=cut

package Bio::EnsEMBL::Exon;

use strict;

use Bio::EnsEMBL::Feature;
use Bio::Seq; # exons have to have sequences...

use Bio::EnsEMBL::Utils::Exception qw( warning throw deprecate );
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Scalar qw( assert_ref );
use Bio::EnsEMBL::DBSQL::SupportingFeatureAdaptor;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Feature);


=head2 new

  Arg [-SLICE]: Bio::EnsEMBL::SLice - Represents the sequence that this
                feature is on. The coordinates of the created feature are
                relative to the start of the slice.
  Arg [-START]: The start coordinate of this feature relative to the start
                of the slice it is sitting on.  Coordinates start at 1 and
                are inclusive.
  Arg [-END]  : The end coordinate of this feature relative to the start of
                the slice it is sitting on.  Coordinates start at 1 and are
                inclusive.
  Arg [-STRAND]: The orientation of this feature.  Valid values are 1,-1,0.
  Arg [-SEQNAME] : (optional) A seqname to be used instead of the default name
                of the of the slice.  Useful for features that do not have an
                attached slice such as protein features.
  Arg [-dbID]   : (optional) internal database id
  Arg [-ADAPTOR]: (optional) Bio::EnsEMBL::DBSQL::BaseAdaptor
  Arg [-PHASE]    : the phase. 
  Arg [-END_PHASE]: the end phase
  Arg [-STABLE_ID]: (optional) the stable id of the exon
  Arg [-VERSION]  : (optional) the version
  Arg [-CREATED_DATE] : (optional) the created date
  Arg [-MODIFIED_DATE]: (optional) the last midifeid date

  Example    : none
  Description: create an Exon object
  Returntype : Bio::EnsEMBL::Exon
  Exceptions : if phase is not valid (i.e. 0,1, 2 -1)
  Caller     : general
  Status     : Stable

=cut

sub new {
  my $class = shift;

  $class = ref $class || $class;

  my $self = $class->SUPER::new( @_ );

  my ( $phase, $end_phase, $stable_id, $version, $created_date,
    $modified_date, $is_current, $is_constitutive )
    = rearrange( [
      "PHASE",        "END_PHASE",
      "STABLE_ID",    "VERSION",
      "CREATED_DATE", "MODIFIED_DATE",
      "IS_CURRENT",   "IS_CONSTITUTIVE"
    ],
    @_
    );

  if ( defined($phase) ) {    # make sure phase is valid.
    $self->phase($phase);
  }

  $self->{'end_phase'}     = $end_phase;
  $self->{'stable_id'}     = $stable_id;
  $self->{'created_date'}  = $created_date;
  $self->{'modified_date'} = $modified_date;

  # Default version
  if ( !defined($version) ) { $version = 1 }
  $self->{'version'} = $version;

  # Default is_current
  if ( !defined($is_current) ) { $is_current = 1 }
  $self->{'is_current'} = $is_current;

  # Default is_constitutive
  if ( !defined($is_constitutive) ) { $is_constitutive = 0 }
  $self->{'is_constitutive'} = $is_constitutive;

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
  Status     : Stable

=cut

sub end_phase {
  my $self = shift;
  if (@_) {
    $self->{'end_phase'} = shift;
  }
  else {
    if ( !defined( $self->{'end_phase'} ) ) {
      my $stable_id = '';
      $stable_id = $self->stable_id if defined $self->stable_id;
      warning("No end phase set in Exon ".$stable_id.". You must set it explicitly.");
    }
  }
  return $self->{'end_phase'};
}


=head2 phase

  Arg [1]    : (optional) int $phase
  Example    :  my $phase = $exon->phase;
                $exon->phase(2);
  Description: Gets/Sets the phase of the exon.
  Returntype : int
  Exceptions : throws if phase is not (0, 1 2 or -1).
  Caller     : general
  Status     : Stable


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
      throw("Bad value ($value) for exon phase. Should only be" .
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
  Status     : Stable

=cut

sub frame {
  my ($self,$value) = @_;

  if( defined $value ) {
    throw("Cannot set frame. Deduced from seq_start and phase");
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

  throw("bad phase in exon ".$self->phase);

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
  Status     : Stable

=cut

sub start {
  my $self = shift;
  # if an arg was provided, flush the internal sequence cache
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
  Status     : Stable

=cut

sub end {
  my $self = shift;
  # if an arg was provided, flush the internal sequence cache
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
  Status     : Stable

=cut

sub strand {
  my $self = shift;
  # if an arg was provided, flush the internal sequence cache
  delete $self->{'_seq_cache'} if(@_);
  return $self->SUPER::strand(@_);
}

=head2 cdna_start

    Arg [1]     : Bio::EnsEMBL::Transcript $transcript
                  The transcript for which cDNA coordinates should be
                  relative to.
    Example     : $cdna_start = $exon->cdna_start($transcript);
    Description : Returns the start position of the exon in cDNA
                  coordinates.
                  Since an exon may be part of one or more transcripts,
                  the relevant transcript must be given as argument to
                  this method.
    Return type : Integer
    Exceptions  : Throws if the given argument is not a transcript.
                  Throws if the first part of the exon maps into a gap.
                  Throws if the exon can not be mapped at all.
    Caller      : General
    Status      : Stable

=cut

sub cdna_start {
  my ($self, $transcript) = @_;
  assert_ref($transcript, 'Bio::EnsEMBL::Transcript', 'transcript');

  my $id = $transcript->dbID();
  
  if(defined $id && exists $self->{cdna_start}->{$id}) {
    return $self->{cdna_start}->{$id};
  }
  
  my $cdna_start;
  my @coords = $transcript->genomic2cdna($self->seq_region_start(), $self->seq_region_end(), $self->strand());
  if(@coords && !$coords[0]->isa('Bio::EnsEMBL::Mapper::Gap')) {
    $cdna_start = $coords[0]->start();
  }
  elsif(@coords) {
    throw "First part of exon maps into gap";
  }
  else {
    throw "Can not map exon";
  }
  
  if(defined $id) {
    $self->{cdna_start}->{$id} = $cdna_start;
  }

  return $cdna_start;
} ## end sub cdna_start

=head2 cdna_end

    Arg [1]     : Bio::EnsEMBL::Transcript $transcript
                  The transcript for which cDNA coordinates should be
                  relative to.
    Example     : $cdna_end = $exon->cdna_end($transcript);
    Description : Returns the end position of the exon in cDNA
                  coordinates.
                  Since an exon may be part of one or more transcripts,
                  the relevant transcript must be given as argument to
                  this method.
    Return type : Integer
    Exceptions  : Throws if the given argument is not a transcript.
                  Throws if the last part of the exon maps into a gap.
                  Throws if the exon can not be mapped at all.
    Caller      : General
    Status      : Stable

=cut

sub cdna_end {
  my ($self, $transcript) = @_;
  assert_ref($transcript, 'Bio::EnsEMBL::Transcript', 'transcript');

  my $id = $transcript->dbID();
  
  if(defined $id && exists $self->{cdna_end}->{$id}) {
    return $self->{cdna_end}->{$id};
  }
  
  my $cdna_end;
  my @coords = $transcript->genomic2cdna($self->seq_region_start(), $self->seq_region_end(), $self->strand());
  if(@coords && !$coords[-1]->isa('Bio::EnsEMBL::Mapper::Gap')) {
    $cdna_end = $coords[-1]->end();
  }
  elsif(@coords) {
    throw "Last part of exon maps into gap";
  }
  else {
    throw "Can not map exon";
  }
  
  if(defined $id) {
    $self->{cdna_end}->{$id} = $cdna_end;
  }

  return $cdna_end;
} ## end sub cdna_end

=head2 cdna_coding_start

    Arg [1]     : Bio::EnsEMBL::Transcript $transcript
                  The transcript for which cDNA coordinates should be
                  relative to.
    Example     : $cdna_coding_start = $exon->cdna_coding_start($transcript);
    Description : Returns the start position of the coding region of the
                  exon in cDNA coordinates.  Returns undef if the whole
                  exon is non-coding.
                  Since an exon may be part of one or more transcripts,
                  the relevant transcript must be given as argument to
                  this method.
    Return type : Integer or undef
    Exceptions  : Throws if the given argument is not a transcript.
    Caller      : General
    Status      : Stable

=cut

sub cdna_coding_start {
  my ($self, $transcript) = @_;
  assert_ref($transcript, 'Bio::EnsEMBL::Transcript', 'transcript');

  my $id = $transcript->dbID();
  
  if(defined $id && exists $self->{cdna_coding_start}->{$id}) {
    return $self->{cdna_coding_start}->{$id};
  }
  
  my $cdna_coding_start;
  my $transcript_coding_start = $transcript->cdna_coding_start();
  if(defined $transcript_coding_start) {
    my $cdna_start = $self->cdna_start($transcript);
    
    if ( $transcript_coding_start < $cdna_start ) {
      # Coding region starts upstream of this exon...

      if ( $transcript->cdna_coding_end() < $cdna_start ) {
        # ... and also ends upstream of this exon.
        $cdna_coding_start = undef;
      }
      else {
        # ... and does not end upstream of this exon.
        $cdna_coding_start = $cdna_start;
      }
    } else {
      # Coding region starts either within or downstream of this
      # exon.

      if ( $transcript_coding_start <= $self->cdna_end($transcript) ) {
        # Coding region starts within this exon.
        $cdna_coding_start = $transcript_coding_start;
      }
      else {
        # Coding region starts downstream of this exon.
        $cdna_coding_start = undef;
      }
    }
  }
  else {
    $cdna_coding_start = undef;
  }
  
  if(defined $id) {
    $self->{cdna_coding_start}->{$id} = $cdna_coding_start;
    $self->{cdna_coding_end}->{$id} = undef if ! defined $cdna_coding_start;
  }
  
  return $cdna_coding_start;
} ## end sub cdna_coding_start

=head2 cdna_coding_end

    Arg [1]     : Bio::EnsEMBL::Transcript $transcript
                  The transcript for which cDNA coordinates should be
                  relative to.
    Example     : $cdna_coding_end = $exon->cdna_coding_end($transcript);
    Description : Returns the end position of the coding region of the
                  exon in cDNA coordinates.  Returns undef if the whole
                  exon is non-coding.
                  Since an exon may be part of one or more transcripts,
                  the relevant transcript must be given as argument to
                  this method.
    Return type : Integer or undef
    Exceptions  : Throws if the given argument is not a transcript.
    Caller      : General
    Status      : Stable

=cut

sub cdna_coding_end {
  my ($self, $transcript) = @_;
  assert_ref($transcript, 'Bio::EnsEMBL::Transcript', 'transcript');

  my $id = $transcript->dbID();
  
  if(defined $id && exists $self->{cdna_coding_end}->{$id}) {
    return $self->{cdna_coding_end}->{$id};
  }
  
  my $cdna_coding_end;
  my $transcript_coding_end = $transcript->cdna_coding_end();
  if(defined $transcript_coding_end) {
    my $cdna_end = $self->cdna_end($transcript);
    
    if ( $transcript_coding_end > $cdna_end ) {
      
      # Coding region ends downstream of this exon...
      if ( $transcript->cdna_coding_start() > $cdna_end ) {
        # ... and also starts downstream of this exon.
        $cdna_coding_end = undef;
      } 
      else {
        # ... and does not start downstream of this exon.
        $cdna_coding_end = $cdna_end;
      }
    }  
    else {
      # Coding region ends either within or upstream of this
      # exon.

      if ( $transcript_coding_end >= $self->cdna_start($transcript) ) {
        # Coding region ends within this exon.
        $cdna_coding_end = $transcript_coding_end;
      } 
      else {
        # Coding region ends upstream of this exon.
        $cdna_coding_end = undef;
      }
    }
  }
  else {
    $cdna_coding_end = undef;
  }
  
  if(defined $id) {
    $self->{cdna_coding_end}->{$id} = $cdna_coding_end;
    $self->{cdna_coding_start}->{$id} = undef if ! defined $cdna_coding_end;
  }
  
  return $cdna_coding_end;
} ## end sub cdna_coding_end

=head2 coding_region_start

    Arg [1]     : Bio::EnsEMBL::Transcript $transcript
    Example     : $coding_region_start =
                    $exon->coding_region_start($transcript);
    Description : Returns the start position of the coding region
                  of the exon in slice-relative coordinates on the
                  forward strand.  Returns undef if the whole exon is
                  non-coding.
                  Since an exon may be part of one or more transcripts,
                  the relevant transcript must be given as argument to
                  this method.
    Return type : Integer or undef
    Exceptions  : Throws if the given argument is not a transcript.
    Caller      : General
    Status      : Stable

=cut

# The implementation of this method is analogous to the implementation
# of cdna_coding_start().

sub coding_region_start {
  my ($self, $transcript) = @_;
  assert_ref($transcript, 'Bio::EnsEMBL::Transcript', 'transcript');

  my $id = $transcript->dbID();
  
  if(defined $id && exists $self->{coding_region_start}->{$id}) {
    return $self->{coding_region_start}->{$id};
  }
  
  my $coding_region_start;
  my $transcript_coding_start = $transcript->coding_region_start();
  if(defined $transcript_coding_start) {
    my $start = $self->start();
    
    if ( $transcript_coding_start < $start ) {
      # Coding region starts upstream of this exon...

      if ( $transcript->coding_region_end() < $start ) {
        # ... and also ends upstream of this exon.
        $coding_region_start = undef;
      } 
      else {
        # ... and does not end upstream of this exon.
        $coding_region_start = $start;
      }
    } 
    else {
      # Coding region starts either within or downstream of this
      # exon.

      if ( $transcript_coding_start <= $self->end() ) {
        # Coding region starts within this exon.
        $coding_region_start = $transcript_coding_start;
      } 
      else {
        # Coding region starts downstream of this exon.
        $coding_region_start = undef;
      }
    }
  }
  else {
    $coding_region_start = undef;
  }
  
  if(defined $id) {
    $self->{coding_region_start}->{$id} = $coding_region_start;
    $self->{coding_region_end}->{$id} = undef if ! defined $coding_region_start;
  }
  
  return $coding_region_start;
} ## end sub coding_region_start

=head2 coding_region_end

    Arg [1]     : Bio::EnsEMBL::Transcript $transcript
    Example     : $coding_region_end =
                    $exon->coding_region_end($transcript);
    Description : Returns the end position of the coding region of
                  the exon in slice-relative coordinates on the
                  forward strand.  Returns undef if the whole exon is
                  non-coding.
                  Since an exon may be part of one or more transcripts,
                  the relevant transcript must be given as argument to
                  this method.
    Return type : Integer or undef
    Exceptions  : Throws if the given argument is not a transcript.
    Caller      : General
    Status      : Stable

=cut

# The implementation of this method is analogous to the implementation
# of cdna_coding_end().

sub coding_region_end {
  my ($self, $transcript) = @_;
  assert_ref($transcript, 'Bio::EnsEMBL::Transcript', 'transcript');

  my $id = $transcript->dbID();
  
  if(defined $id && exists $self->{coding_region_end}->{$id}) {
    return $self->{coding_region_end}->{$id};
  }
  
  my $coding_region_end;
  my $transcript_coding_end = $transcript->coding_region_end();
  if(defined $transcript_coding_end) {
    
    my $end = $self->end();
    if($transcript_coding_end > $end) {
      # Coding region ends downstream of this exon...

      if ( $transcript->coding_region_start() > $end ) {
        # ... and also starts downstream of this exon.
        $coding_region_end = undef;
      } 
      else {
        # ... and does not start downstream of this exon.
        $coding_region_end = $end;
      }
    }
    else {
      # Coding region ends either within or upstream of this
      # exon.
      if ( $transcript_coding_end >= $self->start() ) {
        $coding_region_end = $transcript_coding_end;
      }
      else {
        $coding_region_end = undef;
      }
    }
  }
  else {
    # This is a non-coding transcript.
    $coding_region_end = undef;
  }
  
  if(defined $id) {
    $self->{coding_region_end}->{$id} = $coding_region_end;
    $self->{coding_region_start}->{$id} = undef if ! defined $coding_region_end;
  }
  
  return $coding_region_end;
} ## end sub coding_region_end


=head2 rank

    Arg [1]     : Bio::EnsEMBL::Transcript $transcript
                  The transcript for which the exon rank
                  is requested.
    Example     : $rank = $exon->rank($transcript);
    Description : Returns the rank of the exon relative to
                  the transcript.
                  Since an exon may be part of one or more transcripts,
                  the relevant transcript must be given as argument to
                  this method.
    Return type : Integer
    Exceptions  : Throws if the given argument is not a transcript.
                  Throws if the exon does not belong to the transcript.
    Caller      : General
    Status      : Stable

=cut

sub rank {
  my ($self, $transcript) = @_;
  assert_ref($transcript, 'Bio::EnsEMBL::Transcript', 'transcript');

  my $rank = $transcript->exon_rank($self);

  return $rank;
} ## end sub rank

=head2 slice

  Arg [1]    : Bio::EnsEMBL::Slice
  Example    : $slice = $exon->slice();
  Description: Getter/Setter for the slice this exon is on.  The superclass
               implmentation is overridden to flush the internal sequence
               cache if this value is altered
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub slice {
  my ( $self, $slice ) = @_;

  if ( defined($slice) ) {
    # If a new slice was provided, flush the internal sequence cache and
    # transfer all supporting evidence to the new slice.

    delete $self->{'_seq_cache'};

    if ( exists( $self->{'_supporting_evidence'} ) ) {
      my @new_features;

      for my $old_feature ( @{ $self->{'_supporting_evidence'} } ) {

        my $new_feature;

        if ( defined( $old_feature->slice() ) ) {
          $new_feature = $old_feature->transfer($slice);
        } else {
          # If the old feature does not have a slice, assume transfer is
          # not necessary.
          $new_feature = $old_feature;
        }

        push( @new_features, $new_feature );
      }

      $self->{'_supporting_evidence'} = \@new_features;
    }

    return $self->SUPER::slice($slice);
  } elsif ( @_ > 1 ) {
    return $self->SUPER::slice(undef);
  } else {
    return $self->SUPER::slice();
  }
} ## end sub slice

=head2 equals

  Arg [1]       : Bio::EnsEMBL::Exon exon
  Example       : if ($exonA->equals($exonB)) { ... }
  Description   : Compares two exons for equality.
                  The test for eqality goes through the following list
                  and terminates at the first true match:

                  1. If Bio::EnsEMBL::Feature::equals() returns false,
                     then the exons are *not* equal.
                  2. If both exons have stable IDs: if these are the
                     same, the exons are equal, otherwise not.
                  3. If the exons have the same start, end, strand, and
                     phase, then they are equal, otherwise not.

  Return type   : Boolean (0, 1)

  Exceptions    : Thrown if a non-transcript is passed as the argument.

=cut

sub equals {
  my ( $self, $exon ) = @_;

  if ( !defined($exon) ) { return 0 }
  if ( $self eq $exon ) { return 1 }

  assert_ref( $exon, 'Bio::EnsEMBL::Exon' );

  my $feature_equals = $self->SUPER::equals($exon);
  if ( defined($feature_equals) && $feature_equals == 0 ) {
    return 0;
  }

  if ( defined( $self->stable_id() ) && defined( $exon->stable_id() ) )
  {
    if ( $self->stable_id() eq $exon->stable_id() ) {
      return 1;
    }
    else {
      return 0;
    }
  }

  if ( $self->start() == $exon->start() &&
       $self->end() == $exon->end() &&
       $self->strand() == $exon->strand() &&
       $self->phase() == $exon->phase() &&
       $self->end_phase() == $exon->end_phase() )
  {
    return 1;
  }

  return 0;
} ## end sub equals

=head2 move

  Arg [1]    : int start
  Arg [2]    : int end
  Arg [3]    : (optional) int strand
  Example    : None
  Description: Sets the start, end and strand in one call rather than in 
               3 seperate calls to the start(), end() and strand() methods.
               This is for convenience and for speed when this needs to be
               done within a tight loop.  This overrides the superclass
               move() method so that the internal sequence cache can be
               flushed if the exon if moved.
  Returntype : none
  Exceptions : Thrown is invalid arguments are provided
  Caller     : general
  Status     : Stable

=cut

sub move {
  my $self = shift;
  # flush the internal sequence cache
  delete $self->{'_seq_cache'};
  return $self->SUPER::move(@_);
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
  if( !@_  || ( ref $_[0] && 
         ($_[0]->isa( "Bio::EnsEMBL::Slice" ) or $_[0]->isa( "Bio::EnsEMBL::LRGSlice" ))
        )) {
    deprecate('Calling transform without a coord system name is deprecated.');
    return $self->_deprecated_transform(@_);
  }

  my $new_exon = $self->SUPER::transform( @_ );
  if (not defined $new_exon or
      $new_exon->length != $self->length) {
    return undef;
  }

  if( exists $self->{'_supporting_evidence'} ) {
    my @new_features;
    for my $old_feature ( @{$self->{'_supporting_evidence'}} ) {
      my $new_feature = $old_feature->transform( @_ );
      if (defined $new_feature) {
        push( @new_features, $new_feature );
      }
    }
    $new_exon->{'_supporting_evidence'} = \@new_features;
  }

  #dont want to share the same sequence cache
  delete $new_exon->{'_seq_cache'};

  return $new_exon;
}


=head2 transfer

  Arg [1]    : Bio::EnsEMBL::Slice $destination_slice
  Example    : none
  Description: Moves this Exon to given target slice coordinates. If Features
               are attached they are moved as well. Returns a new exon.
  Returntype : Bio::EnsEMBL::Exon
  Exceptions : none
  Caller     : general
  Status     : Stable

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

  #dont want to share the same sequence cache
  delete $new_exon->{'_seq_cache'};

  return $new_exon;
}


=head2 add_supporting_features

  Arg [1]    : Bio::EnsEMBL::Feature $feature
  Example    : $exon->add_supporting_features(@features);
  Description: Adds a list of supporting features to this exon. 
               Duplicate features are not added.  
               If supporting features are added manually in this
               way, prior to calling get_all_supporting_features then the
               get_all_supporting_features call will not retrieve supporting
               features from the database.
  Returntype : none
  Exceptions : throw if any of the features are not Feature
               throw if any of the features are not in the same coordinate
               system as the exon
  Caller     : general
  Status     : Stable

=cut

sub add_supporting_features {
  my ($self,@features) = @_;

  return unless @features;

  $self->{_supporting_evidence} ||= []; 
  
  # check whether this feature object has been added already
  FEATURE: foreach my $feature (@features) {
    unless($feature && $feature->isa("Bio::EnsEMBL::Feature")) {
      throw("Supporting feat [$feature] not a " .
            "Bio::EnsEMBL::Feature");
    } 
    
    if ((defined $self->slice() && defined $feature->slice())&&
      ( $self->slice()->name() ne $feature->slice()->name())){
      throw("Supporting feat not in same coord system as exon\n" .
            "exon is attached to [".$self->slice()->name()."]\n" .
            "feat is attached to [".$feature->slice()->name()."]");
    }

    foreach my $added_feature ( @{ $self->{_supporting_evidence} } ){
      # compare objects
      if ( $feature == $added_feature ){
  # this feature has already been added
  next FEATURE;
      }
    }
    
    # no duplicate was found, add the feature
    push(@{$self->{_supporting_evidence}},$feature);
  }
}


=head2 flush_supporting_features

  Example     : $exon->flush_supporting_features;
  Description : Removes all supporting evidence from the exon.
  Return type : (Empty) listref
  Exceptions  : none
  Caller      : general
  Status      : Stable

=cut

sub flush_supporting_features {
  my $self = shift;
  $self->{'_supporting_evidence'} = [];
}


=head2 get_all_supporting_features

  Arg [1]    : none
  Example    : @evidence = @{$exon->get_all_supporting_features()};
  Description: Retrieves any supporting features added manually by 
               calls to add_supporting_features. If no features have been
               added manually and this exon is in a database (i.e. it has
               an adaptor), fetch from the database
  Returntype : listreference of Bio::EnsEMBL::BaseAlignFeature objects 
  Exceptions : none
  Caller     : general
  Status     : Stable

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

# This method is only for genebuild backwards compatibility.
# Avoid using it if possible

  Arg [1]    : Bio::EnsEMBL::Feature $features
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
  Status     : Medium Risk

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
      if ($f->entire_seq()->name eq $self->slice()->name) {
  if ($f->end >= $self->start && $f->start <= $self->end && $f->strand == $self->strand) {
    $self->add_supporting_features($f);
  }
      }
    }
  }
}


=head2 stable_id

  Arg [1]    : string $stable_id
  Example    : none
  Description: get/set for attribute stable_id
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub stable_id {
  my $self = shift;
  $self->{'stable_id'} = shift if( @_ );
  return $self->{'stable_id'};
}


=head2 created_date

  Arg [1]    : string $created_date
  Example    : none
  Description: get/set for attribute created_date
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub created_date {
  my $self = shift;
  $self->{'created_date'} = shift if ( @_ );
  return $self->{'created_date'};
}


=head2 modified_date

  Arg [1]    : string $modified_date
  Example    : none
  Description: get/set for attribute modified_date
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub modified_date {
  my $self = shift;
  $self->{'modified_date'} = shift if ( @_ );
  return $self->{'modified_date'};
}


=head2 version

  Arg [1]    : string $version
  Example    : none
  Description: get/set for attribute version
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub version {
   my $self = shift;
  $self->{'version'} = shift if( @_ );
  return $self->{'version'};
}

=head2 stable_id_version

  Arg [1]    : (optional) String - the stable ID with version to set
  Example    : $exon->stable_id("ENSE0000000001.3");
  Description: Getter/setter for stable id with version for this exon.
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub stable_id_version {
    my $self = shift;
    if(my $stable_id = shift) {
	# See if there's an embedded period, assume that's a
	# version, might not work for some species but you
	# should use ->stable_id() and version() if you're worried
	# about ambiguity
	my $vindex = rindex($stable_id, '.');
	# Set the stable_id and version pair depending on if
	# we found a version delimiter in the stable_id
	($self->{stable_id}, $self->{version}) = ($vindex > 0 ?
						  (substr($stable_id,0,$vindex), substr($stable_id,$vindex+1)) :
						  $stable_id, undef);
    }
    return $self->{stable_id} . ($self->{version} ? ".$self->{version}" : '');
}

=head2 is_current

  Arg [1]    : Boolean $is_current
  Example    : $exon->is_current(1)
  Description: Getter/setter for is_current state of this exon.
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub is_current {
  my ( $self, $value ) = @_;

  if ( defined($value) ) {
    $self->{'is_current'} = $value;
  }
  return $self->{'is_current'};
}

=head2 is_constitutive

  Arg [1]    : Boolean $is_constitutive
  Example    : $exon->is_constitutive(0)
  Description: Getter/setter for is_constitutive state of this exon.
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub is_constitutive {
  my ( $self, $value ) = @_;

  if ( defined($value) ) {
    $self->{'is_constitutive'} = $value;
  }
  return $self->{'is_constitutive'};
}

=head2 is_coding

  Arg [1]    : Boolean $is_coding
  Arg [2]    : Bio::EnsEMBL::Transcript
  Example    : $exon->is_coding()
  Description: Says if the exon is within the translation or not
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub is_coding {
  my ( $self, $transcript) = @_;

  if (!$transcript->translate) { return 0; }
  if ($transcript->coding_region_start < $self->start && $self->start < $transcript->coding_region_end) { return 1; }
  if ($transcript->coding_region_end > $self->end && $self->end > $transcript->coding_region_start) { return 1; }

  # to handle cases where transcript coding region can fall within the exon start and exon end, eg: if it is one exon transcript
  if ($transcript->coding_region_start > $self->start &&  $transcript->coding_region_end < $self->end ) { return 1; }
  return 0;
}


=head2 adjust_start_end

  Arg  1     : int $start_adjustment
  Arg  2     : int $end_adjustment
  Example    : none
  Description: returns a new Exon with this much shifted coordinates
  Returntype : Bio::EnsEMBL::Exon
  Exceptions : none
  Caller     : Transcript->get_all_translateable_Exons()
  Status     : Stable

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
  Status     : Stable

=cut

sub peptide {
  my $self = shift;
  my $tr = shift;

  unless($tr && ref($tr) && $tr->isa('Bio::EnsEMBL::Transcript')) {
    throw("transcript arg must be Bio::EnsEMBL:::Transcript not [$tr]");
  }
  
  #convert exons coordinates to peptide coordinates
  my $tmp_exon = $self->transfer($tr->slice);
  if (!$tmp_exon) {
    throw("Couldn't transfer exon to transcript's slice");
  }
  
  my @coords = 
    $tr->genomic2pep($tmp_exon->start, $tmp_exon->end, $tmp_exon->strand);
  
  #filter out gaps
  @coords = grep {$_->isa('Bio::EnsEMBL::Mapper::Coordinate')} @coords;
  
  #if this is UTR then the peptide will be empty string
  my $pep_str = '';
  
  
  if(scalar(@coords) > 1) {    
    my $coord = $self->_merge_ajoining_coords(\@coords);
    if($coord) {
      @coords = ($coord);
    }
    else {
      my ($e_id, $tr_id) = ($self->stable_id(), $tr->stable_id());
      throw("Error. Exon maps to multiple locations in peptide and those". 
      " locations are not continuous." .
      " Is this exon [$e_id] a member of this transcript [$tr_id]?");
    }
  } 
  elsif(scalar(@coords) == 1) {
    my $c = $coords[0];
    my $pep = $tr->translate;
    
    #bioperl doesn't give back residues for incomplete codons
    #make sure we don't subseq too far...
    my ($start, $end);
    $end = ($c->end > $pep->length) ? $pep->length : $c->end; 
    $start = ($c->start < $end) ? $c->start : $end;
    $pep_str = $tr->translate->subseq($start, $end);
  }
  
  return
    Bio::Seq->new( -seq      => $pep_str,
                   -moltype  => 'protein',
                   -alphabet => 'protein',
                   -id       => $self->display_id );
}

=head2 _merge_ajoining_coords
 
   Arg [1]     : ArrayRef of Bio::EnsEMBL::Mapper::Coordinate objects
   Example     : 
   Description : Merges coords which are ajoining or overlapping
   Returntype  : Bio::EnsEMBL::Mapper::Coordinate or undef if it cannot happen
   Exceptions  : Exception if the cooords cannot be condensed into one location
   Caller      : internal
   Status      : Development

=cut
 
sub _merge_ajoining_coords {
  my ($self, $coords) = @_;
  
  my $okay = 1;
  my $coord = shift @{$coords};
  my $start = $coord->start();
  my $last_end = $coord->end();
  foreach my $other_coord (@{$coords}) {
    if( ($last_end + 1) >= $other_coord->start() ) {
      $last_end = $other_coord->end();
    }
    else {
      $okay = 0;
      last;
    }
  }
  
  if(!$okay) {
    return;
  }
  
  my $new_coord = Bio::EnsEMBL::Mapper::Coordinate->new(
  $coord->id(), $start, $last_end, $coord->strand(), $coord->rank());
  return $new_coord;
}




=head2 seq

  Arg [1]    : none
  Example    : my $seq_str = $exon->seq->seq;
  Description: Retrieves the dna sequence of this Exon.
               Returned in a Bio::Seq object.  Note that the sequence may
               include UTRs (or even be entirely UTR).
  Returntype : Bio::Seq or undef
  Exceptions : warning if argument passed,
               warning if exon does not have attatched slice
               warning if exon strand is not defined (or 0)
  Caller     : general
  Status     : Stable

=cut

sub seq {
  my ( $self, $arg ) = @_;

  if ( defined $arg ) {
    warning("seq setting on Exon not supported currently");
    $self->{'_seq_cache'} = $arg->seq();
  }

  if ( !defined( $self->{'_seq_cache'} ) ) {
    my $seq;

    if ( !defined $self->slice() ) {
      warning("Cannot retrieve seq for exon without slice\n");
      return undef;
    }

    if ( !$self->strand() ) {
      warning("Cannot retrieve seq for unstranded exon\n");
      return undef;
    }

    if ($self->slice->is_circular() ) {
  if (    $self->slice->start > $self->slice->end) {
# Normally exons overlapping chromosome origin will have negative feature start, but slice will be from 1 .. length 
# But in case you got an exon attached to  a sub slice try this
      my $mid_point = $self->slice()->seq_region_length() - $self->slice()->start() + 1;
      my $seq1 =  $self->slice()->subseq( $self->start(), $mid_point, $self->strand() );

      my $seq2 = $self->slice()->subseq( $mid_point + 1, $self->end(), $self->strand() );

      $seq = $self->strand() > 0 ? "$seq1$seq2" : "$seq2$seq1";
  } elsif ( $self->start < 0 ||  $self->start > $self->end) {
# Normally exons overlapping chromosome origin will be 0 based, and can have negative start
# But if you go via sub_Slice it gives you chromosome based coordinates, i.e it will have start greater then end
      my $start_point = $self->slice->seq_region_length + $self->slice->start;
      my $mid_point = $self->slice->seq_region_length;
      my $seq1 =  $self->slice->subseq( $self->start, $mid_point, $self->strand);
      my $seq2 = $self->slice->subseq(1, $self->end, $self->strand );
      $seq = $self->strand > 0 ? "$seq1$seq2" : "$seq2$seq1";
  } else {
# End this is the case for genes not overlapping the origin
      $seq = $self->slice()->subseq( $self->start(), $self->end(), $self->strand() );
  }
    } else {
  $seq = $self->slice()->subseq( $self->start(), $self->end(), $self->strand() );
    }

    $self->{'_seq_cache'} = $seq;
  } ## end if ( !defined( $self->...))

  return
    Bio::Seq->new( -seq      => $self->{'_seq_cache'},
                   -id       => $self->display_id,
                   -moltype  => 'dna',
                   -alphabet => 'dna' );
} ## end sub seq


=head2 hashkey

  Arg [1]    : none
  Example    : if(exists $hash{$exon->hashkey}) { do_something(); }
  Description: Returns a unique hashkey that can be used to uniquely identify
               this exon.  Exons are considered to be identical if they share
               the same seq_region, start, end, strand, phase, end_phase.
               Note that this will consider two exons on different slices
               to be different, even if they actually are not. 
  Returntype : string formatted as slice_name-start-end-strand-phase-end_phase
  Exceptions : thrown if not all the necessary attributes needed to generate
               a unique hash value are set
               set
  Caller     : general
  Status     : Stable

=cut

sub hashkey {
  my $self = shift;

  my $slice      = $self->{'slice'}; 
  my $slice_name = ($slice) ? $slice->name() : undef;
  my $start      = $self->{'start'};
  my $end        = $self->{'end'};
  my $strand     = $self->{'strand'};
  my $phase      = $self->{'phase'};
  my $end_phase  = $self->{'end_phase'};

  if(!defined($slice_name)) {
    throw('Slice must be set to generate correct hashkey.');
  }

  if(!defined($start)) {
    warning("start attribute must be defined to generate correct hashkey.");
  }

  if(!defined($end)) {
    throw("end attribute must be defined to generate correct hashkey.");
  }

  if(!defined($strand)) {
    throw("strand attribute must be defined to generate correct hashkey.");
  }

  if(!defined($phase)) {
    throw("phase attribute must be defined to generate correct hashkey.");
  }

  if(!defined($end_phase)) {
    throw("end_phase attribute must be defined to generate correct hashkey.");
  }

  return "$slice_name-$start-$end-$strand-$phase-$end_phase";
}


=head2 display_id

  Arg [1]    : none
  Example    : print $exons->display_id();
  Description: This method returns a string that is considered to be
               the 'display' identifier. For exons this is (depending on
               availability and in this order) the stable Id, the dbID or an
               empty string.
  Returntype : string
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub display_id {
  my $self = shift;
  return $self->{'stable_id'} || $self->dbID || '';
}


=head2 load

  Args          : None
  Example       : $exon->load();
  Description   : The Ensembl API makes extensive use of
                  lazy-loading.  Under some circumstances (e.g.,
                  when copying genes between databases), all data of
                  an object needs to be fully loaded.  This method
                  loads the parts of the object that are usually
                  lazy-loaded.
  Returns       : Nothing.

=cut

sub load {
  my ($self) = @_;

  $self->analysis();
  $self->stable_id();
  $self->get_all_supporting_features();
}

=head2 summary_as_hash

  Example       : $exon_summary = $exon->summary_as_hash();
  Description   : Extends Feature::summary_as_hash
                  Retrieves a summary of this Exon.
  Returns       : hashref of descriptive strings
  Status        : Intended for internal use
=cut

sub summary_as_hash {
  my $self = shift;
  my $summary_ref = $self->SUPER::summary_as_hash;
  $summary_ref->{'constitutive'} = $self->is_constitutive;
  $summary_ref->{'ensembl_phase'} = $self->phase;
  $summary_ref->{'ensembl_end_phase'} = $self->end_phase;
  $summary_ref->{'exon_id'} = $summary_ref->{'id'};
  return $summary_ref;
}


1;
