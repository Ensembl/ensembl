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

Bio::EnsEMBL::Utils::CenteredIntervalTree::Node

=head1 SYNOPSIS


=head1 DESCRIPTION

=head1 METHODS

=cut

package Bio::EnsEMBL::Utils::Tree::Interval::Immutable::Node;

use strict;

use Scalar::Util qw(looks_like_number);
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw);

=head2 new

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($x_center, $s_center, $left, $right) = @_;
  throw 'Node takes arguments (x_center, s_center, left_child, right_child)'
    unless $x_center and $s_center; # and $left and $right;

  throw 'x_center must be a number' unless looks_like_number($x_center);
  assert_ref($s_center, 'ARRAY');
  $left && assert_ref($left, 'Bio::EnsEMBL::Utils::Tree::Interval::Immutable::Node');
  $right && assert_ref($right, 'Bio::EnsEMBL::Utils::Tree::Interval::Immutable::Node');

  my $self = bless({ 'x_center' => $x_center,
		     's_center_beg' => [ map { $_->[1] } sort { $a->[0] <=> $b->[0] } map { [ $_->start, $_ ] } @{$s_center} ],
		     's_center_end' => [ map { $_->[1] } sort { $b->[0] <=> $a->[0] } map { [ $_->end, $_ ] } @{$s_center} ],
		     'lchild'   => $left,
		     'rchild'   => $right }, $class);
  return $self;
}

sub x_center {
  return shift->{x_center};
}

sub s_center_beg {
  return shift->{s_center_beg};
}

sub s_center_end {
  return shift->{s_center_end};
}

sub left {
  return shift->{lchild};
}

sub right {
  return shift->{rchild};
}

1;

