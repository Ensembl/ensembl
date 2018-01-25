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

Bio::EnsEMBL::Utils::Interval

=head1 SYNOPSIS


=head1 DESCRIPTION

=head1 METHODS

=cut

package Bio::EnsEMBL::Utils::Interval;

use strict;

use Scalar::Util qw(looks_like_number);
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw);

=head2 new

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($start, $end, $data) = @_;
  throw 'Must specify interval boundaries [start, end]'
    unless defined $start and defined $end;
  throw 'start must be <= end' if $start > $end;
  
  my $self = bless({ start => $start, end => $end, data => $data }, $class);
  return $self;
}

=head2 start

=cut

sub start {
  my $self = shift;

  return $self->{start};
}

=head2 end

=cut

sub end {
  my $self = shift;

  return $self->{end};
}

=head2 data

=cut

sub data {
  return shift->{data};
}

=head2 is_empty

=cut

sub is_empty {
  my $self = shift;

  return $self->start >= $self->end;
}

=head2 is_point

Determines if the current interval is a single point.

=cut

sub is_point {
  my $self = shift;

  return $self->start == $self->end;
}

=head2 contains

Determines if the current instance contains the query point

=cut

sub contains {
  my ($self, $point) = @_;

  return 0 if $self->is_empty or not defined $point;
  throw 'point must be a number' unless looks_like_number($point);
  
  return ($point >= $self->start and $point <= $self->end);
}

=head2 intersects

=cut

sub intersects {
  my ($self, $interval) = @_;
  assert_ref($interval, 'Bio::EnsEMBL::Utils::Interval');
    
  return ($self->start <= $interval->end and $interval->start <= $self->end);
}

=head2 is_right_of

Checks if this current interval is entirely to the right of a point. More formally,
the method will return true, if for every point x from the current interval the inequality
x > point holds.

=cut

sub is_right_of {
  my ($self, $other) = @_;

  return 0 unless defined $other;
  return $self->start > $other if looks_like_number($other);
  return $self->start > $other->end;
}

=head2 is_left_of

Checks if this current interval is entirely to the left of another interval. More formally,
the method will return true, if for every point x from the current interval the inequality
x < point holds.

=cut

sub is_left_of {
  my ($self, $other) = @_;

  return 0 unless defined $other;
  return $self->end < $other if looks_like_number($other);
  return $self->end < $other->start;
}

1;

