=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2017] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Utils::CenteredIntervalTree

=head1 SYNOPSIS


=head1 DESCRIPTION

Heavily inspired by https://github.com/tylerkahn/intervaltree-python

=head1 METHODS

=cut

package Bio::EnsEMBL::Utils::CenteredIntervalTree;

use strict;

use Data::Dumper;

use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Utils::CenteredIntervalTree::Node;

=head2 new

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my $intervals = shift;
  $intervals && assert_ref($intervals, 'ARRAY');
  
  my $self = bless({}, $class);

  $self->{top_node} = $self->_divide_intervals($intervals);
  return $self;
}

sub root {
  return shift->{top_node};
}

# sub search {
#   my ($self, $start, $end) = @_;
#   defined $start or throw 'Undefined point or interval';

#   # interval search
#   if (defined $end) {
#     my $result = [];
#     for my $i ($start .. $end + 1) {
#       push @{$result}, @{$self->search($i)};
#     }

#     return sort_by_begin(uniq($result));
#   } 

#   # point search 
#   return $self->_search($self->{top_node}, $start, []);
# }

# sub _search {
#   my ($self, $node, $point, $result) = @_;

#   # can be optimized, see https://en.wikipedia.org/wiki/Interval_tree
#   for my $i (@{$node->scenter}) {
#     push @{$result}, $i if $i->start <= $point and $i->end >= $point;
#   }

#   if ($point < $node->x_center and $node->left_child) {
#     my $left_results = $self->_search($node->left_child, $point, []);
#     push @{$result}, @{$left_results} if $left_results;
#   }

#   if ($point > $node->x_center and $node->right_child) {
#     my $right_results = $self->_search($node->right_child, $point, []);
#     push @{$result}, @{$right_results} if $right_results;
#   } 

#   return uniq($result);
# }

sub query {
  my ($self, $interval) = @_;

  return [] unless $interval;
  return $self->_query_point($self->root, $interval->start, []) if $interval->is_point;

  my $result = [];
  return $result unless $self->root or $interval->is_empty;

  my $node = $self->root;
  while ($node) {
    if ($interval->contains($node->x_center)) {
      push @{$result}, @{$node->s_center_beg};
      $self->_range_query_left($node->left, $interval, $result);
      $self->_range_query_right($node->right, $interval, $result);
      last;
    }
    if ($interval->is_left_of($node->x_center)) {
      foreach my $s_beg (@{$node->s_center_beg}) {
	last unless $interval->intersects($s_beg);
	push @{$result}, $s_beg;
      }
      $node = $node->left;
    } else {
      foreach my $s_end (@{$node->s_center_end}) {
	last unless $interval->intersects($s_end);
	push @{$result}, $s_end;
      }
      $node = $node->right;
    }
  }

  return sort_by_begin(uniq($result));
}

sub _query_point {
  my ($self, $node, $point, $result) = @_;

  return $result unless $node;

  # if x is less than x_center, the leftmost set of intervals, S_left, is considered
  if ($point <= $node->x_center) {
    # From Wikipedia
    # if x is less than x_center, we know that all intervals in S_center end after x,
    # or they could not also overlap x_center. Therefore, we need only find those intervals
    # in S_center that begin before x. We can consult the lists of S_center that have already
    # been constructed. Since we only care about the interval beginnings in this scenario,
    # we can consult the list sorted by beginnings.
    # Suppose we find the closest number no greater than x in this list. All ranges from the
    # beginning of the list to that found point overlap x because they begin before x and end
    # after x (as we know because they overlap x_center which is larger than x).
    # Thus, we can simply start enumerating intervals in the list until the startpoint value exceeds x.
    foreach my $s_beg (@{$node->s_center_beg}) {
      last if $s_beg->is_right_of($point);
      push @{$result}, $s_beg;
    }

    # since x < x_center, we also consider the leftmost set of intervals
    return $self->_query_point($node->left, $point, $result);
  } else {
    # if x is greater than x_center, we know that all intervals in S_center must begin before x,
    # so we find those intervals that end after x using the list sorted by interval endings.
    foreach my $s_end (@{$node->s_center_end}) {
      last if $s_end->is_left_of($point);
      push @{$result}, $s_end;
    }
    
    # since x > x_center, we also consider the rightmost set of intervals
    return $self->_query_point($node->right, $point, $result);
  }
  
  return sort_by_begin(uniq($result));
}

# This corresponds to the left branch of the range search, once we find a node, whose
# midpoint is contained in the query interval. All intervals in the left subtree of that node
# are guaranteed to intersect with the query, if they have an endpoint greater or equal than
# the start of the query interval. Basically, this means that every time we branch to the left
# in the binary search, we need to add the whole right subtree to the result set.

sub _range_query_left {
  my ($self, $node, $interval, $result) = @_;
  
  while ($node) {
    if ($interval->contains($node->x_center)) {
      push @{$result}, @{$node->s_center_beg};
      if ($node->right) {
	# in-order traversal of the right subtree to add all its intervals
	$self->_in_order_traversal($node->right, $result);
      }
      $node = $node->left;
    } else {
      foreach my $seg_end (@{$node->s_center_end}) {
	last if $seg_end->is_left_of($interval);
	push @{$result}, $seg_end;
      }
      $node = $node->right;
    }
  }
}

# This corresponds to the right branch of the range search, once we find a node, whose
# midpoint is contained in the query interval. All intervals in the right subtree of that node
# are guaranteed to intersect with the query, if they have an endpoint smaller or equal than
# the end of the query interval. Basically, this means that every time we branch to the right
# in the binary search, we need to add the whole left subtree to the result set.

sub _range_query_right {
  my ($self, $node, $interval, $result) = @_;

  while ($node) {
    if ($interval->contains($node->x_center)) {
      push @{$result}, @{$node->s_center_beg};
      if ($node->left) {
	# in-order traversal of the left subtree to add all its intervals
	$self->_in_order_traversal($node->left, $result);
      }
      $node = $node->right;
    } else {
      foreach my $seg_beg (@{$node->s_center_beg}) {
	last if $seg_beg->is_right_of($interval);
	push @{$result}, $seg_beg;
      }
      $node = $node->left;
    }
  }
}

sub in_order_traversal {
  my ($self) = @_;

  my $result = [];
  $self->_in_order_traversal($self->root, $result);

  return $result;
}

sub _in_order_traversal {
  my ($self, $node, $result) = @_;

  return unless $node;
  $result ||= [];

  $self->_in_order_traversal($node->left, $result);
  push @{$result}, @{$node->s_center_beg};
  $self->_in_order_traversal($node->right, $result);
}

sub _divide_intervals {
  my ($self, $intervals) = @_;

  return undef unless scalar @{$intervals};

  my $x_center = $self->_center($intervals);
  my ($s_center, $s_left, $s_right) = ([], [], []);
  
  foreach my $interval (@{$intervals}) {
    if ($interval->end < $x_center) {
      push @{$s_left}, $interval;
    } elsif ($interval->start > $x_center) {
      push @{$s_right}, $interval;
    } else {
      push @{$s_center}, $interval;
    }
  }

  my $node = Bio::EnsEMBL::Utils::CenteredIntervalTree::Node->new($x_center,
								  $s_center,
								  $self->_divide_intervals($s_left),
								  $self->_divide_intervals($s_right));
}

sub _center {
  my ($self, $intervals) = @_;

  my $sorted_intervals = sort_by_begin($intervals);
  my $len = scalar @{$sorted_intervals};
  
  return $sorted_intervals->[int($len/2)]->start;
}

sub sort_by_begin {
  my $intervals = shift;

  return [ map { $_->[1] } sort { $a->[0] <=> $b->[0] } map { [ $_->start, $_ ] } @{$intervals} ];
}

sub uniq {
  my $intervals = shift;

  use Tie::RefHash;
  tie my %seen, 'Tie::RefHash';
  
  return [ grep { ! $seen{ $_ }++ } @{$intervals} ];
}

1;

