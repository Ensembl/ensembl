package Bio::EnsEMBL::IdMapping::SyntenyRegion;

=head1 NAME


=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS


=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http://www.ensembl.org/info/about/code_licence.html for details.

=head1 AUTHOR

Patrick Meidl <meidl@ebi.ac.uk>, Ensembl core API team

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<ensembl-dev@ebi.ac.uk>

=cut


use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::Utils::Exception qw(throw warning);


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


sub source_start {
  my $self = shift;
  $self->[0] = shift if (@_);
  return $self->[0];
}


sub source_end {
  my $self = shift;
  $self->[1] = shift if (@_);
  return $self->[1];
}


sub source_strand {
  my $self = shift;
  $self->[2] = shift if (@_);
  return $self->[2];
}


sub source_seq_region_name {
  my $self = shift;
  $self->[3] = shift if (@_);
  return $self->[3];
}


sub target_start {
  my $self = shift;
  $self->[4] = shift if (@_);
  return $self->[4];
}


sub target_end {
  my $self = shift;
  $self->[5] = shift if (@_);
  return $self->[5];
}


sub target_strand {
  my $self = shift;
  $self->[6] = shift if (@_);
  return $self->[6];
}


sub target_seq_region_name {
  my $self = shift;
  $self->[7] = shift if (@_);
  return $self->[7];
}


sub score {
  my $self = shift;
  $self->[8] = shift if (@_);
  return $self->[8];
}


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

#
# extend this SyntenyRegion to span a $factor * $score more area
#
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

