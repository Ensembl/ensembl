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

use strict;
use warnings;

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
