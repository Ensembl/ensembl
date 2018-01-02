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

Bio::EnsEMBL::FeaturePair - Stores sequence Features which are
themselves hits to other sequence features.

=head1 SYNOPSIS

    my $feat = Bio::EnsEMBL::FeaturePair->new(
      -start      => 132_231,
      -end        => 132_321,
      -strand     => -1,
      -slice      => $slice,
      -hstart     => 10,
      -hend       => 100,
      -hstrand    => 1,
      -score      => 100,
      -percent_id => 92.0,
      -hseqname   => 'ALUSX10.1',
      -analysis   => $analysis
    );

    my $hit_start  = $feat->hstart();
    my $hit_end    = $feat->hend();
    my $hit_strand = $feat->hstrand();
    my $analysis   = $feat->analysis();

=head1 DESCRIPTION

A sequence feature object where the feature is itself a feature on
another sequence - e.g. a blast hit where residues 1-40 of a protein
sequence SW:HBA_HUMAN has hit to bases 100 - 220 on a genomic sequence
HS120G22.  The genomic sequence coordinates are represented by the
start, end, strand attributes while the protein (hit) coordinates are
represented by the hstart, hend, hstrand attributes.

  $clone = $slice_adpator->fetch_by_region( 'clone', 'HS120G22' );

  $fp = Bio::EnsEMBL::FeaturePair(
    -start      => 100,
    -end        => 220,
    -strand     => 1,
    -slice      => $clone,
    -hstart     => 1,
    -hend       => 40,
    -hstrand    => 1,
    -percent_id => 92.0,
    -score      => 100,
    -hseqname   => 'SW:HBA_HUMAN',
    -species    => 'Homo sapiens',
    -hspecies   => 'Homo sapiens'
  );

=head1 METHODS

=cut

package Bio::EnsEMBL::FeaturePair;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate warning);

@ISA = qw(Bio::EnsEMBL::Feature);

=head2 new

  Arg [HSTART]    : int - The start of the hit region (optional)
  Arg [HEND]      : int - The end of the hit region (optional)
  Arg [HSTRAND]   : (0,1,-1) - The strand of the hit region (optional)
  Arg [PERCENT_ID]: float -  The precentage identity of the hit (optional)
  Arg [SCORE]     : float -  The score of the hit (optional)
  Arg [HSEQNAME]  : string - The name of the hit sequence (optional)
  Arg [P_VALUE]   : float -  The pvalue or evalue (optional)
  Arg [SPECIES]   : string - The species the query sequence is from (optional)
  Arg [HSPECIES]  : string - The species the hit sequence is from (optional)
  Arg [COVERAGE]  : string - The % of the query that this feature pair covers
  Arg [HCOVERAGE] : string - The % of the target this this feature pair covers
  Arg [EXTRA_DATA]: HashRef - Additional data, specified as name, value attribute pairs (optional)
  Arg [...]       : Named superclass constructor args (Bio::EnsEMBL::Feature)
  Example    : $feat = Bio::EnsEMBL::FeaturePair->new(-start    => 132_231,
                                              -end      => 132_321,
                                              -strand   => -1,
                                              -slice    => $slice,
                                              -hstart   => 10,
                                              -hend     => 100,
                                              -hstrand  => 1,
                                              -score    => 100,
                                              -percent_id => 92.0,
                                              -hseqname => 'ALUSX10.1',
                                              -analysis => $analysis);
  Description: Creates a new Bio::EnsEMBL::FeaturePair object
  Returntype : Bio::EnsEMBL::FeaturePair
  Exceptions : throw if start > end
               throw if invalid strand is provided
  Caller     : general
  Status     : Stable

=cut

sub new {
  my $caller = shift;

  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  my ($hstart, $hend, $hstrand, $percent_id, $score, $species, $hspecies, $p_value, $hseqname, $f1, $f2, $coverage, $hcoverage, $group_id, $level_id, $external_db_id, $extra_data, $external_db_name, $external_display_db_name, $hdescription) = rearrange(['HSTART', 'HEND', 'HSTRAND', 'PERCENT_ID', 'SCORE', 'SPECIES', 'HSPECIES', 'P_VALUE', 'HSEQNAME', 'FEATURE1', 'FEATURE2', 'COVERAGE', 'HCOVERAGE', 'GROUP_ID', 'LEVEL_ID', 'EXTERNAL_DB_ID', 'EXTRA_DATA', 'DBNAME', 'DB_DISPLAY_NAME', 'HDESCRIPTION'], @_);

  if (defined($hstart) && defined($hend) && ($hend < $hstart)) {
	throw('HSTART must be less than or equal to HEND');
  }

  if (defined($hstrand) && $hstrand != 1 && $hstrand != -1 && $hstrand != 0) {
	throw('HSTRAND must be one of (0,1,-1)');
  }

  $self->{'hstart'}          = $hstart;
  $self->{'hend'}            = $hend;
  $self->{'hstrand'}         = $hstrand;
  $self->{'score'}           = $score;
  $self->{'percent_id'}      = $percent_id;
  $self->{'species'}         = $species;
  $self->{'hspecies'}        = $hspecies;
  $self->{'hseqname'}        = $hseqname;
  $self->{'coverage'}        = $coverage;
  $self->{'hcoverage'}       = $hcoverage;
  $self->{'p_value'}         = $p_value;
  $self->{'group_id'}        = $group_id;
  $self->{'level_id'}        = $level_id;
  $self->{'external_db_id'}  = $external_db_id;
  $self->{'extra_data'}      = $extra_data;
  $self->{'dbname'}          = $external_db_name;
  $self->{'db_display_name'} = $external_display_db_name;
  $self->{'hdescription'}    = $hdescription;
  #
  # Feature1 and Feature2 arg handling for backwards compatibility
  #
  if ($f1) {
	deprecate("Using FEATURE1 arg to construct FeaturePairs" . " is deprecated.\nUse the args START,END,STRAND,SLICE instead");

	#eval because we are not exactly sure what f1 arg will look like
	eval {
	  $self->{'start'}    = $f1->start();
	  $self->{'end'}      = $f1->end();
	  $self->{'strand'}   = $f1->strand();
	  $self->{'slice'}    = $f1->contig();
	  $self->{'analysis'} = $f1->analysis() if ($f1->analysis());
	};
  }

  if ($f2) {
	deprecate("Using FEATURE2 arg to construct FeaturePairs is deprecated" . "\nUse the args HSTART,HEND,HSTRAND,HSEQNAME instead");

	#eval because we are not exactly sure what f2 arg will look like
	eval {
	  $self->{'hseqname'} = $f2->seqname();
	  $self->{'hstart'}   = $f2->start();
	  $self->{'hend'}     = $f2->end();
	  $self->{'hstrand'}  = $f2->strand();
	  $self->{'analysis'} = $f2->analysis() if ($f2->analysis());
	};
  }

  return $self;
} ## end sub new

=head2 hseqname

  Arg [1]    : string $hseqname (optional)
  Example    : $hseqname = $fp->hseqname();
  Description: Getter/Setter for the name of the hit sequence
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub hseqname {
  my $self = shift;
  $self->{'hseqname'} = shift if (@_);
  return $self->{hseqname};
}

=head2 hstart

  Arg [1]    : string $hstart (optional)
  Example    : $hstart = $fp->hstart();
  Description: Getter/Setter for the start coordinate on the hit sequence
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub hstart {
  my $self = shift;
  $self->{'hstart'} = shift if (@_);
  return $self->{'hstart'};
}

=head2 hend

  Arg [1]    : string $hend (optional)
  Example    : $hend = $fp->hend();
  Description: Getter/Setter for the end coordinate on the hit sequence
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub hend {
  my $self = shift;
  $self->{'hend'} = shift if (@_);
  return $self->{'hend'};
}

=head2 hstrand

  Arg [1]    : int $hstrand (optional)
  Example    : $hstrand = $fp->hstrand
  Description: Getter/Setter for the orientation of the hit on the hit sequence
  Returntype : 0,1,-1
  Exceptions : thrown 
  Caller     : general
  Status     : Stable

=cut

sub hstrand {
  my $self = shift;

  if (@_) {
	my $hstrand = shift;
	if (defined($hstrand) && $hstrand != 1 && $hstrand != 0 && $hstrand != -1) {
	  throw('hstrand must be one of (-1,0,1)');
	}
	$self->{'hstrand'} = $hstrand;
  }

  return $self->{'hstrand'};
}

=head2 hslice

  Arg [1]    : (optional) Bio::EnsEMBL::Slice $slice
  Example    : $hseqname = $featurepair->hslice()->seq_region_name();
  Description: Getter/Setter for the Slice that is associated with this 
               hit feature.  The slice represents the underlying sequence that this
               feature is on.  Note that this method call is analagous to the
               old SeqFeature methods contig(), entire_seq(), attach_seq(),
               etc.
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : thrown if an invalid argument is passed
  Caller     : general
  Status     : Stable

=cut

sub hslice {
  my $self = shift;

  if (@_) {
	my $sl = shift;
	if (defined($sl) && (!ref($sl) || !($sl->isa('Bio::EnsEMBL::Slice')))) {
	  throw('slice argument must be a Bio::EnsEMBL::Slice');
	}

	$self->{'hslice'} = $sl;
  }

  return $self->{'hslice'};
}

=head2 hseq_region_name

  Arg [1]    : none
  Example    : print $feature->hseq_region_name();
  Description: Gets the name of the hseq_region which this feature is on.
               Returns undef if this Feature is not on a hslice.
  Returntype : string or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub hseq_region_name {
  my $self  = shift;
  my $slice = $self->{'hslice'};

  return ($slice) ? $slice->seq_region_name() : undef;
}

=head2 hseq_region_strand

  Arg [1]    : none
  Example    : print $feature->hseq_region_strand();
  Description: Returns the strand of the hseq_region which this feature is on 
               (i.e. feature_strand * slice_strand)
               Returns undef if this Feature is not on a hslice.
  Returntype : 1,0,-1 or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub hseq_region_strand {
  my $self  = shift;
  my $slice = $self->{'hslice'};

  return ($slice) ? $slice->strand()*$self->{'hstrand'} : undef;
}

=head2 hseq_region_start

  Arg [1]    : none
  Example    : print $feature->hseq_region_start();
  Description: Convenience method which returns the absolute start of this
               feature on the hseq_region, as opposed to the relative (hslice) 
               position.

               Returns undef if this feature is not on a hslice.
  Returntype : int or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub hseq_region_start {
  my $self  = shift;
  my $slice = $self->{'hslice'};

  return undef if (!$slice);

  if ($slice->strand == 1) {
	return undef if (!defined($self->{'hstart'}));
	return $slice->start() + $self->{'hstart'} - 1;
  } else {
	return undef if (!defined($self->{'hend'}));
	return $slice->end() - $self->{'hend'} + 1;
  }
}

=head2 hseq_region_end

  Arg [1]    : none
  Example    : print $feature->hseq_region_end();
  Description: Convenience method which returns the absolute end of this
               feature on the hseq_region, as opposed to the relative (hslice)
               position.

               Returns undef if this feature is not on a hslice.
  Returntype : int or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub hseq_region_end {
  my $self  = shift;
  my $slice = $self->{'hslice'};

  return undef if (!$slice);

  if ($slice->strand == 1) {
	return undef if (!defined($self->{'hend'}));
	return $slice->start() + $self->{'hend'} - 1;
  } else {
	return undef if (!defined($self->{'hstart'}));
	return $slice->end() - $self->{'hstart'} + 1;
  }
}

=head2 score

  Arg [1]    : float $score (optional)
  Example    : $score = $fp->score();
  Description: Getter/Setter for the score of this feature pair
  Returntype : float
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub score {
  my $self = shift;
  $self->{'score'} = shift if (@_);
  return $self->{'score'};
}

=head2 percent_id

  Arg [1]    : float $percent_id (optional)
  Example    : $percent_id = $fp->percent_id();
  Description: Getter/Setter for the percentage identity of this feature pair
  Returntype : float
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub percent_id {
  my $self = shift;
  $self->{'percent_id'} = shift if (@_);
  return $self->{'percent_id'};
}

=head2 species

 Arg [1]    : string $genus_species_name (optional)
              e.g. Homo_sapiens or Mus_musculus
 Example    : $species = $fp->species();
 Description: get/set on the species of feature1
 Returntype : string
 Execeptions: none
 Caller     : general
 Status     : Stable

=cut

sub species {
  my $self = shift;
  $self->{'species'} = shift if (@_);
  return $self->{'species'};
}

=head2 hspecies

 Arg [1]    : string $genus_species_name (optional)
              e.g. Homo_sapiens or Mus_musculus
 Example    : $hspecies = $fp->hspecies
 Description: get/set on the species of feature2
 Returntype : string
 Execeptions: none
 Caller     : general
 Status     : Stable

=cut

sub hspecies {
  my $self = shift;
  $self->{'hspecies'} = shift if (@_);
  return $self->{'hspecies'};
}

=head2 coverage

  Arg [1]    : number (percentage) $coverage (optional)
  Example    : $cov = $fp->coverage();
  Description: Getter/Setter for the % of the query covered by the feature
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub coverage {
  my $self = shift;
  $self->{'coverage'} = shift if (@_);
  return $self->{'coverage'};
}

=head2 hcoverage

  Arg [1]    : number (percentage) $hcoverage (optional)
  Example    : $hcov = $fp->hcoverage();
  Description: Getter/Setter for the % of the target covered by the feature
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub hcoverage {
  my $self = shift;
  $self->{'hcoverage'} = shift if (@_);
  return $self->{'hcoverage'};
}

=head2 external_db_id

  Arg [1]    : int  $external_db_id (optional)
  Example    : $ex_db = $fp->external_db_id();
  Description: Getter/Setter for the external_db_id taregt source database feature
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub external_db_id {
  my $self = shift;
  $self->{'external_db_id'} = shift if (@_);
  return $self->{'external_db_id'};
}

=head2 db_name

  Arg [1]    : string  $external_db_name (optional)
  Example    : $ex_db_name = $fp->dbname();
  Description: Getter/Setter for the external_db_name attribute, name of external database
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub db_name {
  my $self = shift;
  $self->{'dbname'} = shift if (@_);
  return $self->{'dbname'};
}

=head2 db_display_name

  Arg [1]    : string  $db_display_name (optional)
  Example    : $ex_db_display_name = $fp->db_display_name();
  Description: Getter/Setter for the db_display_name attribute 
               The preferred display name for the external database. 
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub db_display_name {
  my $self = shift;
  $self->{'db_display_name'} = shift if (@_);
  return $self->{'db_display_name'};
}

=head2 p_value

  Arg [1]    : float $p_value (optional)
  Example    : $eval = $fp->p_value
  Description: Getter Setter for the evalue / pvalue of this feature
  Returntype : float
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub p_value {
  my $self = shift;
  $self->{'p_value'} = shift if (@_);
  return $self->{'p_value'};
}

=head2 hdescription

  Arg [1]    : String (optional)
  Example    : $des = $fp->hdescription()
  Description: Getter Setter for optional description of this feature
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub hdescription {
  my $self = shift;
  $self->{'hdescription'} = shift if (@_);
  return $self->{'hdescription'};
}

=head2 display_id

  Arg [1]    : none
  Example    : print $fp->display_id();
  Description: This method returns a string that is considered to be
               the 'display' identifier.  For feature pairs this is the 
               hseqname if it is available otherwise it is an empty string.
  Returntype : string
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub display_id {
  my $self = shift;
  return $self->{'hseqname'} || '';
}

=head2 identical_matches

 Arg [1]    : int $identical_matches (optional)
 Example    : 
 Description: get/set on the number of identical matches
 Returntype : int
 Execeptions: none
 Caller     : general
  Status     : Stable

=cut

sub identical_matches {
  my ($self, $arg) = @_;

  if (defined($arg)) {
	return $self->{'_identical_matches'} = $arg;
  }
  return $self->{'_identical_matches'};
}

=head2 positive_matches

 Arg [1]    : int $positive_matches (optional)
 Example    : 
 Description: get/set on the number of positive matches
 Returntype : int
 Execeptions: none
 Caller     : general
  Status     : Stable

=cut

sub positive_matches {
  my ($self, $arg) = @_;

  if (defined($arg)) {
	return $self->{'_positive_matches'} = $arg;
  }
  return $self->{'_positive_matches'};
}

=head2 group_id
 
  Arg [1]    : int $group_id
  Example    : none
  Description: get/set for attribute group_id
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable
 
=cut

sub group_id {
  my ($self, $arg) = @_;

  if (defined $arg) {
	$self->{'group_id'} = $arg;
  }
  return $self->{'group_id'};
}

=head2 level_id
 
  Arg [1]    : int $level_id
  Example    : none
  Description: get/set for attribute level_id
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable
 
=cut

sub level_id {
  my ($self, $arg) = @_;

  if (defined $arg) {
	$self->{'level_id'} = $arg;
  }
  return $self->{'level_id'};
}

=head1 DEPRECATED METHODS


=head2 invert

  Arg [1]    : (optional) Bio::EnsEMBL::Slice $newslice
  Example    : $feature->invert();
  Description: This method is used to swap the hit and query sides of this
               feature in place.  A new slice may optionally provided which
               this feature will be placed on.  If no slice is provided the
               feature slice will be set to undef.
  Returntype : none
  Exceptions : none
  Caller     : pipeline (BlastMiniGenewise)

=cut

sub invert {
  my ($self, $slice) = @_;

  if (!defined $slice && defined $self->hslice) {
	$slice = $self->hslice;
  }

  my $hstart       = $self->{'hstart'};
  my $hend         = $self->{'hend'};
  my $hstrand      = $self->{'hstrand'};
  my $hspecies     = $self->{'hspecies'};
  my $hseqname     = $self->{'hseqname'};
  my $hdescription = $self->{'hdescription'};

  my $start   = $self->{'start'};
  my $end     = $self->{'end'};
  my $strand  = $self->{'strand'};
  my $species = $self->{'species'};
  my $seqname = $self->seqname();

  $self->{'start'}   = $hstart;
  $self->{'end'}     = $hend;
  $self->{'strand'}  = $hstrand;
  $self->{'species'} = $hspecies;
  $self->{'seqname'} = $hseqname if (defined($hseqname));

  $self->{'hstart'}   = $start;
  $self->{'hend'}     = $end;
  $self->{'hstrand'}  = $strand;
  $self->{'hseqname'} = $seqname;
  $self->{'hspecies'} = $species;

  $self->{'hdescription'} = $hdescription;

  $self->{'hslice'} = $self->slice;
  $self->{'slice'}  = $slice;
} ## end sub invert


sub extra_data {
  my $self = shift;
  $self->{'extra_data'} = shift if (@_);
  return $self->{'extra_data'};
}

sub type {
  my $self = shift;
  $self->{'extra_data'}->{'type'} = shift if (@_);
  if (exists $self->{'extra_data'}) {
	return $self->{'extra_data'}->{'type'};
  }
  return;
}

1;
