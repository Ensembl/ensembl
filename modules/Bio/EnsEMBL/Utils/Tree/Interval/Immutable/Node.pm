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

Bio::EnsEMBL::Utils::Tree::Interval::Immutable::Node

=head1 DESCRIPTION

Instances of this class represent nodes in a centered immutable interval tree.

=head1 METHODS

=cut

package Bio::EnsEMBL::Utils::Tree::Interval::Immutable::Node;

use strict;

use Scalar::Util qw(looks_like_number);
use Bio::EnsEMBL::Utils::Exception qw(throw);

=head2 new

  Arg [1]     : scalar, $x_center
                The mid point of the intervals represented by the node
  Arg [2]     : Arrayref of Bio::EnsEMBL::Utils::Interval instances
                Represents the set of intervals overlapping $x_center
  Arg [3]     : Bio::EnsEMBL::Utils::Tree::Interval::Immutable::Node
                The left subtree
  Arg [4]     : Bio::EnsEMBL::Utils::Tree::Interval::Immutable::Node
                The right subtree
  Description : Constructor. Creates a new centered immutable interval tree node instance
  Returntype  : Bio::EnsEMBL::Utils::Tree::Interval::Immutable::Node
  Exceptions  : if $x_center is not a number or the remaining args are not of the
                correct type
  Caller      : general

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($x_center, $s_center, $left, $right) = @_;
  throw 'Node takes arguments (x_center, s_center, left_child, right_child)'
    unless $x_center and $s_center; # and $left and $right;

  throw 'x_center must be a number' unless looks_like_number($x_center);

  # Bio::EnsEMBL::Utils::Scalar::assert_ref is high-overhead, and was
  # the source of a significant performance penalty in pipelines where
  # millions of new nodes are being instantiated. Replaced with simpler
  # low-overhead type validation.
  throw 's_center must be an array' unless (ref($s_center) eq 'ARRAY');

  if(defined($left)) {
      unless (ref($left) eq 'Bio::EnsEMBL::Utils::Tree::Interval::Immutable::Node') {
          throw 'left must be a Bio::EnsEMBL::Utils::Tree::Interval::Immutable::Node';
      }
  }

  if(defined($right)) {
      unless (ref($right) eq 'Bio::EnsEMBL::Utils::Tree::Interval::Immutable::Node') {
          throw 'right must be a Bio::EnsEMBL::Utils::Tree::Interval::Immutable::Node';
      }
  }

  my $self = bless({ 'x_center' => $x_center,
		     's_center_beg' => [ map { $_->[1] } sort { $a->[0] <=> $b->[0] } map { [ $_->start, $_ ] } @{$s_center} ],
		     's_center_end' => [ map { $_->[1] } sort { $b->[0] <=> $a->[0] } map { [ $_->end, $_ ] } @{$s_center} ],
		     'lchild'   => $left,
		     'rchild'   => $right }, $class);
  return $self;
}

=head2 x_center

  Arg []      : none
  Description : Return the center point of the intervals represented by the node
  Returntype  : scalar
  Exceptions  : none
  Caller      : general

=cut

sub x_center {
  return shift->{x_center};
}

=head2 s_center_beg

  Arg []      : none
  Description : Returns the set of intervals containing the center point sorted by their start point
  Returntype  : Arrayref of Bio::EnsEMBL::Utils::Interval
  Exceptions  : none
  Caller      : general

=cut

sub s_center_beg {
  return shift->{s_center_beg};
}

=head2 s_center_end

  Arg []      : none
  Description : Returns the set of intervals containing the center point sorted by their end point
  Returntype  : Arrayref of Bio::EnsEMBL::Utils::Interval
  Exceptions  : none
  Caller      : general

=cut

sub s_center_end {
  return shift->{s_center_end};
}

=head2 left

  Arg []      : none
  Description : Returns the node's left child
  Returntype  : Bio::EnsEMBL::Utils::Tree::Interval::Immutable::Node
  Exceptions  : none
  Caller      : general

=cut

sub left {
  return shift->{lchild};
}

=head2 right

  Arg []      : none
  Description : Returns the node's right child
  Returntype  : Bio::EnsEMBL::Utils::Tree::Interval::Immutable::Node
  Exceptions  : none
  Caller      : general

=cut

sub right {
  return shift->{rchild};
}

1;

