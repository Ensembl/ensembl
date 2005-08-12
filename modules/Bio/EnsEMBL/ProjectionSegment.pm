#
# Ensembl module for Bio::EnsEMBL::ProjectionSegment
#
#
# Copyright Team Ensembl
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::ProjectionSegment - part of the list that is returned from
project function calls

=head1 SYNOPSIS


   $slice = $sa->fetch_by_region('chromosome', 'X', 1_000_000, 2_000_000);

   my $projection = $slice->project( "clone" );

   foreach my $projection_segment ( @$projection ) {

     print "  from_start ",$projection_segment->from_start(), "\n";
     print "  from_end   ",$projection_segment->from_end(), "\n";
     print "  to_Slice   ",$projection_segment->to_Slice()->name(), "\n";
   }


=head1 DESCRIPTION

The ProjectionSegment is a helper object to make the arrays returned by 
project more accessible. Instead of writing $segment->[0], $segment->[1] or
$segment->[2] its possible to use the more descriptive notation of
$segment->from_start(), $segement->from_end(), $segment->to_Slice().

=head1 AUTHOR - Arne Stabenau

=head1 CONTACT

This modules is part of the Ensembl project http://www.ensembl.org

Questions can be posted to the ensembl-dev mailing list:
ensembl-dev@ebi.ac.uk

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
