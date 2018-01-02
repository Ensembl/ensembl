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

Bio::EnsEMBL::IdMapping::SyntenyRegion - object representing syntenic regions

=head1 SYNOPSIS

  # create a new SyntenyRegion from a source and a target gene
  my $sr = Bio::EnsEMBL::IdMapping::SyntenyRegion->new_fast( [
      $source_gene->start,  $source_gene->end,
      $source_gene->strand, $source_gene->seq_region_name,
      $target_gene->start,  $target_gene->end,
      $target_gene->strand, $target_gene->seq_region_name,
      $entry->score,
  ] );

  # merge with another SyntenyRegion
  my $merged_sr = $sr->merge($sr1);

  # score a gene pair against this SyntenyRegion
  my $score =
    $sr->score_location_relationship( $source_gene1, $target_gene1 );

=head1 DESCRIPTION

This object represents a synteny between a source and a target location.
SyntenyRegions are built from mapped genes, and the their score is
defined as the score of the gene mapping. For merged SyntenyRegions,
scores are combined.

=head1 METHODS

  new_fast
  source_start
  source_end
  source_strand
  source_seq_region_name
  target_start
  target_end
  target_strand
  target_seq_region_name
  score
  merge
  stretch
  score_location_relationship
  to_string

=cut

package Bio::EnsEMBL::IdMapping::SyntenyRegion;

use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::Utils::Exception qw(throw warning);


=head2 new_fast

  Arg[1]      : Arrayref $array_ref - the arrayref to bless into the
                SyntenyRegion object 
  Example     : my $sr = Bio::EnsEMBL::IdMapping::SyntenyRegion->new_fast([
                  ]);
  Description : Constructor. On instantiation, source and target regions are
                reverse complemented so that source is always on forward strand.
  Return type : a Bio::EnsEMBL::IdMapping::SyntenyRegion object
  Exceptions  : none
  Caller      : Bio::EnsEMBL::IdMapping::SyntenyFramework
  Status      : At Risk
              : under development

=cut

sub new_fast {
  my $class = shift;
  my $array_ref = shift;

  # reverse complement source and target so that source is always on forward
  # strand; this will make merging and other comparison operations easier
  # at later stages
  if ($array_ref->[2] == -1) {
    $array_ref->[2] = 1;
    $array_ref->[6] = -1 * $array_ref->[6];
  }
  
  return bless $array_ref, $class;
}


=head2 source_start

  Arg[1]      : (optional) Int - source location start coordinate
  Description : Getter/setter for source location start coordinate.
  Return type : Int
  Exceptions  : none
  Caller      : Bio::EnsEMBL::IdMapping::SyntenyFramework
  Status      : At Risk
              : under development

=cut

sub source_start {
  my $self = shift;
  $self->[0] = shift if (@_);
  return $self->[0];
}


=head2 source_end

  Arg[1]      : (optional) Int - source location end coordinate
  Description : Getter/setter for source location end coordinate.
  Return type : Int
  Exceptions  : none
  Caller      : Bio::EnsEMBL::IdMapping::SyntenyFramework
  Status      : At Risk
              : under development

=cut

sub source_end {
  my $self = shift;
  $self->[1] = shift if (@_);
  return $self->[1];
}


=head2 source_strand

  Arg[1]      : (optional) Int - source location strand
  Description : Getter/setter for source location strand.
  Return type : Int
  Exceptions  : none
  Caller      : Bio::EnsEMBL::IdMapping::SyntenyFramework
  Status      : At Risk
              : under development

=cut

sub source_strand {
  my $self = shift;
  $self->[2] = shift if (@_);
  return $self->[2];
}


=head2 source_seq_region_name

  Arg[1]      : (optional) String - source location seq_region name
  Description : Getter/setter for source location seq_region name.
  Return type : String
  Exceptions  : none
  Caller      : Bio::EnsEMBL::IdMapping::SyntenyFramework
  Status      : At Risk
              : under development

=cut

sub source_seq_region_name {
  my $self = shift;
  $self->[3] = shift if (@_);
  return $self->[3];
}


=head2 target_start

  Arg[1]      : (optional) Int - target location start coordinate
  Description : Getter/setter for target location start coordinate.
  Return type : Int
  Exceptions  : none
  Caller      : Bio::EnsEMBL::IdMapping::SyntenyFramework
  Status      : At Risk
              : under development

=cut

sub target_start {
  my $self = shift;
  $self->[4] = shift if (@_);
  return $self->[4];
}


=head2 target_end

  Arg[1]      : (optional) Int - target location end coordinate
  Description : Getter/setter for target location end coordinate.
  Return type : Int
  Exceptions  : none
  Caller      : Bio::EnsEMBL::IdMapping::SyntenyFramework
  Status      : At Risk
              : under development

=cut

sub target_end {
  my $self = shift;
  $self->[5] = shift if (@_);
  return $self->[5];
}


=head2 target_strand

  Arg[1]      : (optional) Int - target location strand
  Description : Getter/setter for target location strand.
  Return type : Int
  Exceptions  : none
  Caller      : Bio::EnsEMBL::IdMapping::SyntenyFramework
  Status      : At Risk
              : under development

=cut

sub target_strand {
  my $self = shift;
  $self->[6] = shift if (@_);
  return $self->[6];
}


=head2 target_seq_region_name

  Arg[1]      : (optional) String - target location seq_region name
  Description : Getter/setter for target location seq_region name.
  Return type : String
  Exceptions  : none
  Caller      : Bio::EnsEMBL::IdMapping::SyntenyFramework
  Status      : At Risk
              : under development

=cut

sub target_seq_region_name {
  my $self = shift;
  $self->[7] = shift if (@_);
  return $self->[7];
}


=head2 score

  Arg[1]      : (optional) Float - score
  Description : Getter/setter for the score between source and target location.
  Return type : Int
  Exceptions  : none
  Caller      : Bio::EnsEMBL::IdMapping::SyntenyFramework
  Status      : At Risk
              : under development

=cut

sub score {
  my $self = shift;
  $self->[8] = shift if (@_);
  return $self->[8];
}


=head2 merge

  Arg[1]      : Bio::EnsEMBL::IdMapping::SyntenyRegion $sr - another
                SyntenyRegion
  Example     : $merged_sr = $sr->merge($other_sr);
  Description : Merges two overlapping SyntenyRegions if they meet certain
                criteria (see documentation in the code for details). Score is
                calculated as a combined distance score. If the two
                SyntenyRegions aren't mergeable, this method returns undef.
  Return type : Bio::EnsEMBL::IdMapping::SyntenyRegion or undef
  Exceptions  : warns on bad scores
  Caller      : Bio::EnsEMBL::IdMapping::SyntenyFramework
  Status      : At Risk
              : under development

=cut

sub merge {
  my ($self, $sr) = @_;

  # must be on same seq_region
  if ($self->source_seq_region_name ne $sr->source_seq_region_name or
      $self->target_seq_region_name ne $sr->target_seq_region_name) {
    return 0;
  }

  # target must be on same strand
  return 0 unless ($self->target_strand == $sr->target_strand);

  # find the distance of source and target pair and compare
  my $source_dist = $sr->source_start - $self->source_start;
  my $target_dist;
  if ($self->target_strand == 1) {
    $target_dist = $sr->target_start - $self->target_start;
  } else {
    $target_dist = $self->target_end - $sr->target_end;
  }

  # prevent division by zero error
  if ($source_dist == 0 or $target_dist == 0) {
    warn("WARNING: source_dist ($source_dist) and/or target_dist ($target_dist) is zero.\n");
    return 0;
  }

  # calculate a distance score
  my $dist = $source_dist - $target_dist;
  $dist = -$dist if ($dist < 0);
  my $d1 = $dist/$source_dist;
  $d1 = -$d1 if ($d1 < 0);
  my $d2 = $dist/$target_dist;
  $d2 = -$d2 if ($d2 < 0);
  my $dist_score = 1 - $d1 - $d2;

  # distance score must be more than 50%
  return 0 if ($dist_score < 0.5);

  my $new_score = $dist_score * ($sr->score + $self->score)/2;

  if ($new_score > 1) {
    warn("WARNING: Bad merge score: $new_score\n");
  }

  # extend SyntenyRegion to cover both sources and targets, set merged score
  # and return
  if ($sr->source_start < $self->source_start) {
    $self->source_start($sr->source_start);
  }
  if ($sr->source_end > $self->source_end) {
    $self->source_end($sr->source_end);
  }
  
  if ($sr->target_start < $self->target_start) {
    $self->target_start($sr->target_start);
  }
  if ($sr->target_end > $self->target_end) {
    $self->target_end($sr->target_end);
  }

  $self->score($new_score);

  return $self;
}


=head2 stretch

  Arg[1]      : Float $factor - stretching factor
  Example     : $stretched_sr = $sr->stretch(2);
  Description : Extends this SyntenyRegion to span a $factor * $score more area.
  Return type : Bio::EnsEMBL::IdMapping::SyntenyRegion
  Exceptions  : none
  Caller      : Bio::EnsEMBL::IdMapping::SyntenyFramework
  Status      : At Risk
              : under development

=cut

sub stretch {
  my ($self, $factor) = @_;

  my $source_adjust = int(($self->source_end - $self->source_start + 1) *
    $factor * $self->score);
  $self->source_start($self->source_start - $source_adjust);
  $self->source_end($self->source_end + $source_adjust);
  #warn sprintf("  sss %d %d %d\n", $source_adjust, $self->source_start,
  #  $self->source_end);
  
  my $target_adjust = int(($self->target_end - $self->target_start + 1) *
    $factor * $self->score);
  $self->target_start($self->target_start - $target_adjust);
  $self->target_end($self->target_end + $target_adjust);

  return $self;
}


=head2 score_location_relationship

  Arg[1]      : Bio::EnsEMBL::IdMapping::TinyGene $source_gene - source gene
  Arg[2]      : Bio::EnsEMBL::IdMapping::TinyGene $target_gene - target gene
  Example     : my $score = $sr->score_location_relationship($source_gene,
                  $target_gene);
  Description : This function calculates how well the given source location
                interpolates on given target location inside this SyntenyRegion.

                Scoring is done the following way: Source and target location
                are normalized with respect to this Regions source and target.
                Source range will then be somewhere close to 0.0-1.0 and target
                range anything around that.

                The extend of the covered area between source and target range
                is a measurement of how well they agree (smaller extend is
                better). The extend (actually 2*extend) is reduced by the size
                of the regions. This will result in 0.0 if they overlap
                perfectly and bigger values if they dont.

                This is substracted from 1.0 to give the score. The score is
                likely to be below zero, but is cut off at 0.0f.

                Finally, the score is multiplied with the score of the synteny
                itself.
  Return type : Float
  Exceptions  : warns if score out of range
  Caller      : Bio::EnsEMBL::IdMapping::SyntenyFramework
  Status      : At Risk
              : under development

=cut



sub score_location_relationship {
  my ($self, $source_gene, $target_gene) = @_;

  # must be on same seq_region
  if (($self->source_seq_region_name ne $source_gene->seq_region_name) or
      ($self->target_seq_region_name ne $target_gene->seq_region_name)) {
    return 0;
  }

  # strand relationship must be the same (use logical XOR to find out)
  if (($self->source_strand == $source_gene->strand) xor
      ($self->target_strand == $target_gene->strand)) {
    return 0;
  }

  # normalise source location
  my $source_rel_start = ($source_gene->start - $self->source_start) /
                         ($self->source_end - $self->source_start + 1);

  my $source_rel_end = ($source_gene->end - $self->source_start + 1) /
                       ($self->source_end - $self->source_start + 1);

  #warn "  aaa ".$self->to_string."\n";
  #warn sprintf("  bbb %.6f %.6f\n", $source_rel_start, $source_rel_end);

  # cut off if the source location is completely outside
  return 0 if ($source_rel_start > 1.1 or $source_rel_end < -0.1);
  
  # normalise target location
  my ($target_rel_start, $target_rel_end);
  my $t_length = $self->target_end - $self->target_start + 1;

  if ($self->target_strand == 1) {

    $target_rel_start = ($target_gene->start - $self->target_start) / $t_length;

    $target_rel_end = ($target_gene->end - $self->target_start + 1) / $t_length;

  } else {
    $target_rel_start = ($self->target_end - $target_gene->end) / $t_length;
    $target_rel_end = ($self->target_end - $target_gene->start + 1) / $t_length;
  }

  my $added_range = (($target_rel_end > $source_rel_end) ? $target_rel_end :
                                                           $source_rel_end) -
              (($target_rel_start < $source_rel_start) ? $target_rel_start :
                                                         $source_rel_start);

  my $score = $self->score * (1 - (2 * $added_range - $target_rel_end -
    $source_rel_end + $target_rel_start + $source_rel_start));

  #warn "  ccc ".sprintf("%.6f:%.6f:%.6f:%.6f:%.6f\n", $added_range,
  #  $source_rel_start, $source_rel_end, $target_rel_start, $target_rel_end);

  $score = 0 if ($score < 0);

  # sanity check
  if ($score > 1) {
    warn "Out of range score ($score) for ".$source_gene->id.":".
      $target_gene->id."\n";
  }

  return $score;
}


=head2 to_string

  Example     : print LOG $sr->to_string, "\n";
  Description : Returns a string representation of the SyntenyRegion object.
                Useful for debugging and logging.
  Return type : String
  Exceptions  : none
  Caller      : Bio::EnsEMBL::IdMapping::SyntenyFramework
  Status      : At Risk
              : under development

=cut

sub to_string {
  my $self = shift;
  return sprintf("%s:%s-%s:%s %s:%s-%s:%s %.6f",
    $self->source_seq_region_name,
    $self->source_start,
    $self->source_end,
    $self->source_strand,
    $self->target_seq_region_name,
    $self->target_start,
    $self->target_end,
    $self->target_strand,
    $self->score
  );
}


1;

