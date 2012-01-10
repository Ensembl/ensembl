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

Bio::EnsEMBL::ProjectionSegment - part of the list that is returned from
project function calls

=head1 SYNOPSIS

  $slice =
    $sa->fetch_by_region( 'chromosome', 'X', 1_000_000, 2_000_000 );

  my $projection = $slice->project("clone");

  foreach my $projection_segment (@$projection) {
    print( "  from_start ", $projection_segment->from_start(), "\n" );
    print( "  from_end   ", $projection_segment->from_end(),   "\n" );
    print( "  to_Slice   ",
      $projection_segment->to_Slice()->name(), "\n" );
  }

=head1 DESCRIPTION

The ProjectionSegment is a helper object to make the arrays returned by
project more accessible. Instead of writing $segment->[0], $segment->[1]
or $segment->[2] its possible to use the more descriptive notation of
$segment->from_start(), $segement->from_end(), $segment->to_Slice().

=head1 METHODS

=cut

package Bio::EnsEMBL::ProjectionSegment;

#
# WARNING: THIS CLASS IS REPRESENTED BY A BLESSED ARRAY REFERENCE
#  NOT A HASH REFERENCE
#





=head2 from_start

  Args       : none
  Example    : $coord_in_fetaure_start = $segment->from_start()
  Description: First element in projects returned segment lists
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub from_start {
  my $self = shift;
  return $self->[0];
}



=head2 from_end

  Args       : none
  Example    : $coord_in_feature_end = $segment->from_end()
  Description: Second element in projects returned segment lists
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub from_end {
  my $self = shift;
  return $self->[1];
}




=head2 to_Slice

  Args       : none
  Example    : $target_slice = $segment->to_Slice()
  Description: Third element in projects returned segment lists
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub to_Slice {
  my $self = shift;
  return $self->[2];
}



1;
