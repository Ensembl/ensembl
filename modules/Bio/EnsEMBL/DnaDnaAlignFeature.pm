package Bio::EnsEMBL::DnaDnaAlignFeature;

# EnsEMBL module for storing dna-dna pairwise alignments
#
# Cared for by Michele Clamp <michele@sanger.ac.uk>
#
# You may distribute this module under the same terms as perl itself
#

=head1 NAME

  Bio::EnsEMBL::DnaDnaAlignFeature - Ensembl specific dna-dna pairwise alignment feature

=head1 SYNOPSIS

  See BaseAlignFeature

=cut


use Bio::EnsEMBL::BaseAlignFeature;


use vars qw(@ISA);
use strict;

@ISA = qw( Bio::EnsEMBL::BaseAlignFeature );


=head2 _hit_unit

  Arg [1]    : none
  Example    : none
  Description: PRIVATE implementation of abstract superclass method.  Returns
               1 as the 'unit' used for the hit sequence. 
  Returntype : int
  Exceptions : none
  Caller     : Bio::EnsEMBL::BaseAlignFeature

=cut

sub _hit_unit {
  return 1;
}



=head2 _query_unit

  Arg [1]    : none
  Example    : none
  Description: PRIVATE implementation of abstract superclass method Returns
               1 as the 'unit' used for the hit sequence.
  Returntype : int
  Exceptions : none
  Caller     : Bio::EnsEMBL::BaseAlignFeature

=cut

sub _query_unit {
  return 1;
}

=head2 restrict_between_positions

  Arg  1     : int $start
  Arg  2     : int $end
  Arg  3     : string $sequence_name
               The Feature can either be restricted on Ensembl coords
               or on hit coords. The sequence_name decides which one.
  Example    : none
  Description: Makes a new AlignFeature that sits in the new coords
  Returntype : DnaDnaAlignFeature
  Exceptions : returns undef if it cant do the job
  Caller     : general

=cut

sub restrict_between_positions {
  my ($self,$start,$end,$seqref) = @_;

  unless (defined $start && $start =~ /^\d+$/) {
    $self->throw("The first argument is not defined or is not an integer");
  }
  unless (defined $end && $end =~ /^\d+$/) {
    $self->throw("The second argument is not defined or is not an integer");
  }
  unless (defined $seqref &&
          ($seqref eq "seqname" || $seqref eq "hseqname")) {
    $self->throw("The third argument is not defined or is not equal to 'seqname' or 'hseqname'");
  }

# symbolic method references should be forbidden!

  my ($start_method1,$end_method1,$strand_method1,$start_method2,$end_method2,$strand_method2) =
    qw(start end strand hstart hend hstrand);

  if ($seqref eq "hseqname") {
    ($start_method1,$end_method1,$strand_method1,$start_method2,$end_method2,$strand_method2) =
    qw(hstart hend hstrand start end strand);
  }

  my @restricted_features;

  foreach my $ungapped_feature ($self->ungapped_features) {

    if ($ungapped_feature->$start_method1() > $end ||
        $ungapped_feature->$end_method1() < $start) {

      next;

    } elsif ($ungapped_feature->$end_method1() <= $end &&
             $ungapped_feature->$start_method1() >= $start) {

      push @restricted_features, $ungapped_feature;

    } else {

      if ($ungapped_feature->$strand_method1() eq $ungapped_feature->$strand_method2()) {

        if ($ungapped_feature->$start_method1() < $start) {

          my $offset = $start - $ungapped_feature->$start_method1();
          $ungapped_feature->$start_method1($start);
          $ungapped_feature->$start_method2($ungapped_feature->$start_method2() + $offset);

        }
        if ($ungapped_feature->$end_method1() > $end) {

          my $offset = $ungapped_feature->$end_method1() - $end;
          $ungapped_feature->$end_method1($end);
          $ungapped_feature->$end_method2($ungapped_feature->$end_method2() - $offset);

        }
      } else {

        if ($ungapped_feature->$start_method1() < $start) {

          my $offset = $start - $ungapped_feature->$start_method1();
          $ungapped_feature->$start_method1($start);
          $ungapped_feature->$end_method2($ungapped_feature->$end_method2() - $offset);

        }
        if ($ungapped_feature->$end_method1() > $end) {

          my $offset = $ungapped_feature->$end_method1() - $end;
          $ungapped_feature->$end_method1($end);
          $ungapped_feature->$start_method2($ungapped_feature->$start_method2() + $offset);

        }
      }
      
      push @restricted_features, $ungapped_feature;
    }
  }

  if (scalar @restricted_features) {
    my $DnaDnaAlignFeature = new Bio::EnsEMBL::DnaDnaAlignFeature('-features' =>\@restricted_features);
    if (defined $self->slice) {
      $DnaDnaAlignFeature->slice($self->slice);
    }
    if (defined $self->hslice) {
      $DnaDnaAlignFeature->hslice($self->hslice);
    }
    return $DnaDnaAlignFeature;
  } else {
    return undef;
  }
}

1;
