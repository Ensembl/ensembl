=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Utils::Interval

=head1 SYNOPSIS

  # let's get an interval spanning 9e5 bp and associated it with some data
  my $i1 = Bio::EnsEMBL::Utils::Interval->new(1e5, 1e6, { 'key1' => 'value1', 'key2' => 'value2' });

  # and another one which overlaps with the previous,
  # but with scalar associated data
  my $i2 = Bio::EnsEMBL::Utils::Interval->new(2e5, 3e5, 'a string' );

  warn "Empty interval(s)\n" if $i1->is_empty or $i2->is_empty;
  warn "Point interval(s)\n" if $i1->is_point or $i2->is_point;

  if ($i1->intersects($i2)) {
    print "I1 and I2 overlap\n";
  } else {
    print "I1 and I2 do not overlap\n";
  }

  # If an interval is defined with a start > end, then it is assumed
  # to be spanning the origin on a circular chromosome
  my $i3 = Bio::EnsEMBL::Utilities::Interval->new(1e5, 1e2);
  warn "Interval spans the origin" if $i3->spans_origin;

  etc.

=head1 DESCRIPTION

A class representing an interval defined on a genomic region. Instances of this
class can store arbitrarily defined data. If created with start > end, then it
is assumed that this interval is on a circular chromosome spanning the origin.

=head1 METHODS

=cut

package Bio::EnsEMBL::Utils::Interval;

use strict;

use Scalar::Util qw(looks_like_number);
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw);

=head2 new

  Arg [1]     : scalar $start
                The start coordinate of the region
  Arg [2]     : scalar $end
                The end coordinate of the region
  Arg [3]     : (optional) $data
                The data associated with the interval, can be anything
  Example     : my $i = Bio::EnsEMBL::Utils::Interval(1e2, 2e2, { 'key' => 'value' });
                my $i2 = Bio::EnsEMBL::Utilities::Interval(1e5, 1e2);
                $i->spans_origin # returns 0
                $i2->spans_origin # returns 1
  Description : Constructor. Creates a new instance
  Returntype  : Bio::EnsEMBL::Utils::Interval
  Exceptions  : Throws an exception if start and end are not defined.
  Caller      : general

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($start, $end, $data) = @_;
  throw 'Must specify interval boundaries [start, end]'
    unless defined $start and defined $end;

  my $spans_origin = 0;
  if ($start > $end) {
      $spans_origin = 1;
  }
  
  my $self = bless({ start => $start,
                     end => $end,
                     data => $data ,
                     spans_origin => $spans_origin},
                   $class);
  return $self;
}

=head2 start

  Arg []      : none
  Description : Returns the start coordinate of the region
  Returntype  : scalar
  Exceptions  : none
  Caller      : general

=cut

sub start {
  my $self = shift;

  return $self->{start};
}

=head2 end

  Arg []      : none
  Description : Returns the end coordinate of the region
  Returntype  : scalar
  Exceptions  : none
  Caller      : general

=cut

sub end {
  my $self = shift;

  return $self->{end};
}

=head2 data

  Arg []      : none
  Description : Returns the data associated with the region
  Returntype  : Any
  Exceptions  : none
  Caller      : general

=cut

sub data {
  my $self = shift;

  return $self->{data};
}

=head2 spans_origin

  Arg []      : none
  Description : Returns whether this interval was created spanning zero
                (more particularly: if the interval was instantiated with start > end)
  Returntype  : boolean
  Exceptions  : none
  Caller      : general

=cut

sub spans_origin {
  my $self = shift;

  return $self->{spans_origin};
}

=head2 is_empty

  Arg []      : none
  Description : Returns whether or not the interval is empty
  Returntype  : boolean
  Exceptions  : none
  Caller      : general

=cut

sub is_empty {
  my $self = shift;

  if ($self->spans_origin) {
    return ($self->end >= $self->start);
  } else {
    return ($self->start >= $self->end);
  }
}

=head2 is_point

  Arg []      : none
  Description : Determines if the current interval is a single point
  Returntype  : boolean
  Exceptions  : none
  Caller      : general

=cut

sub is_point {
  my $self = shift;

  return $self->start == $self->end;
}

=head2 contains

  Arg [1]     : scalar, the point coordinate 
  Description : Determines if the current instance contains the query point
  Returntype  : boolean
  Exceptions  : none
  Caller      : general

=cut

sub contains {
  my ($self, $point) = @_;

  return 0 if $self->is_empty or not defined $point;
  throw 'point must be a number' unless looks_like_number($point);
  
  if ($self->spans_origin) {
    return ($point >= $self->start or $point <= $self->end);
  } else {
    return ($point >= $self->start and $point <= $self->end);
  }
}

=head2 intersects

  Arg [1]     : An instance of Bio::EnsEMBL::Utils::Interval
  Description : Determines if the the instance intersects the given interval
  Returntype  : boolean
  Exceptions  : none
  Caller      : general

=cut

sub intersects {
  my ($self, $interval) = @_;
  assert_ref($interval, 'Bio::EnsEMBL::Utils::Interval');

  if ($self->spans_origin and $interval->spans_origin) {
    return 1;
  } elsif ($self->spans_origin or $interval->spans_origin) {
    return ($interval->end >= $self->start or $interval->start <= $self->end);
  } else {
    return ($self->start <= $interval->end and $interval->start <= $self->end);
  }
}

=head2 is_right_of

  Arg [1]     : An instance of Bio::EnsEMBL::Utils::Interval or a scalar
  Description : Checks if this current interval is entirely to the right of a point
                or Interval.
                More formally, the method will return true, if for every point x from 
                the current interval the inequality x > point holds, where point
                is either a single scalar, or point is the end of another Interval.
                If spans_origin is true for either this Interval or an Interval
                passed in, then this method returns false.
  Returntype  : boolean
  Exceptions  : none
  Caller      : general

=cut

sub is_right_of {
  my ($self, $other) = @_;

  return 0 unless defined $other;

  if ( looks_like_number($other) ) {
    return $self->spans_origin ?
           throw "is_right_of not defined for an interval that spans the origin" :
           $self->start > $other;
  } elsif ($self->spans_origin or $other->spans_origin) {
    throw "is_right_of not defined for an interval that spans the origin";
  } else {
    return $self->start > $other->end;
  }
}

=head2 is_left_of

  Arg [1]     : An instance of Bio::EnsEMBL::Utils::Interval or a scalar
  Description : Checks if this current interval is entirely to the left of a point
                or Interval.
                More formally, the method will return true, if for every point x from 
                the current interval the inequality x < point holds, where point
                is either a single scalar, or point is the start of another Interval.
                If spans_origin is true for either this Interval or an Interval
                passed in, then this method returns false
  Returntype  : boolean
  Exceptions  : none
  Caller      : general

=cut

sub is_left_of {
  my ($self, $other) = @_;

  return 0 unless defined $other;

  if ( looks_like_number($other) ) {
    return $self->spans_origin ?
           throw "is_left_of not defined for an interval that spans the origin" :
           $self->end < $other;
  } elsif ($self->spans_origin or $other->spans_origin) {
    throw "is_left_of not defined for an interval that spans the origin";
  } else {
    return $self->end < $other->start;
  }
}

1;

