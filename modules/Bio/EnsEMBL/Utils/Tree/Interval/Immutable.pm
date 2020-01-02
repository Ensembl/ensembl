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

Bio::EnsEMBL::Utils::Tree::Interval::Immutable

=head1 SYNOPSIS

  # define a set of intervals to be added to the tree
  my $intervals = [ Bio::EnsEMBL::Utils::Interval->new(121626874, 122092717),
		    Bio::EnsEMBL::Utils::Interval->new(121637917, 121658918),
		    Bio::EnsEMBL::Utils::Interval->new(122096077, 124088369) ];

  # initialise the tree with the above intervals
  my $tree = Bio::EnsEMBL::Utils::Tree::Interval::Immutable->new($intervals);

  # point query
  my $results = $tree->query(121779004);
  if (scalar @$results) {
    print "Intervals contain 121779004\n";
  }

  # same query, but use interval query
  my $results = $tree->query(121779004, 121779004);
  if (scalar @$results) {
    print "Found containing interval: [", $result->[0]->start, ', ', $result->[0]->end, "\n";
  }

=head1 DESCRIPTION

An implementation of an immutable interval tree. Immutable means the tree is
initialised with a fixed set of intervals at creation time. Intervals cannot
be added to or removed from the tree during its life cycle.

Implementation heavily inspired by https://github.com/tylerkahn/intervaltree-python

This implementation does not support Intervals having a start > end - i.e.
intervals spanning the origin of a circular chromosome.

=head1 METHODS

=cut

package Bio::EnsEMBL::Utils::Tree::Interval::Immutable;

use strict;

use Tie::RefHash;

use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Utils::Tree::Interval::Immutable::Node;
use Bio::EnsEMBL::Utils::Interval;


=head2 new

  Arg [1]     : Arrayref of Bio::EnsEMBL::Utils::Interval instances
  Example     : my $tree = Bio::EnsEMBL::Utils::Tree::Immutable([ $i1, $i2, $i3 ]);
  Description : Constructor. Creates a new immutable tree instance
  Returntype  : Bio::EnsEMBL::Utils::Tree::Interval::Immutable
  Exceptions  : none
  Caller      : general

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my $intervals = shift;
  if (defined $intervals ) {
    assert_ref($intervals, 'ARRAY');
  }
  
  my $self = bless({}, $class);

  $self->{top_node} = $self->_divide_intervals($intervals);

  return $self;
}

=head2 root

  Arg []      : none
  Example     : my $root = $tree->root();
  Description : Return the tree top node
  Returntype  : Bio::EnsEMBL::Utils::Tree::Interval::Immutable::Node
  Exceptions  : none
  Caller      : general

=cut

sub root {
  return shift->{top_node};
}

=head2 query

  Arg [1]     : scalar, $start
                Where the query interval begins
  Arg [2]     : (optional) scalar, $end
                Where the query interval ends
  Example     : my $results = $tree->query(121626874, 122092717);
  Description : Query the tree if its intervals overlap the interval whose start
                and end points are specified by the argument list.
                If end is not specified, it is assumed to be the same as start
                so effectively making a point query.
  Returntype  : An arrayref of Bio::EnsEMBL::Utils::Interval instances
  Exceptions  : none
  Caller      : general

=cut

sub query {
  my ($self, $start, $end) = @_;
  
  my $interval;
  if (defined $start) {
    $end = $start unless defined $end;
    $interval = Bio::EnsEMBL::Utils::Interval->new($start, $end);
  }
  
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
  my ($self, $intervals, $sorted) = @_;

  return undef unless scalar @{$intervals};

  my $sorted_intervals;
  if ($sorted) {
      $sorted_intervals = $intervals;
  } else {
      $sorted_intervals = sort_by_begin($intervals);
  }

  my $x_center = $self->_center_sorted($sorted_intervals);
  my ($s_center, $s_left, $s_right) = ([], [], []);
  
  foreach my $sorted_interval (@{$sorted_intervals}) {
    if ($sorted_interval->spans_origin) {
      throw "Cannot build a tree containing an interval that spans the origin";
    }
    if ($sorted_interval->end < $x_center) {
      push @{$s_left}, $sorted_interval;
    } elsif ($sorted_interval->start > $x_center) {
      push @{$s_right}, $sorted_interval;
    } else {
      push @{$s_center}, $sorted_interval;
    }
  }

  my $node = Bio::EnsEMBL::Utils::Tree::Interval::Immutable::Node->new($x_center,
								       $s_center,
								       $self->_divide_intervals($s_left, 1),
								       $self->_divide_intervals($s_right, 1));
}

sub _center {
  my ($self, $intervals) = @_;

  my $sorted_intervals = sort_by_begin($intervals);
  my $len = scalar @{$sorted_intervals};
  
  return $sorted_intervals->[int($len/2)]->start;
}

sub _center_sorted {
  my ($self, $sorted_intervals) = @_;
  my $len = scalar @{$sorted_intervals};

  return $sorted_intervals->[int($len/2)]->start;
}

sub sort_by_begin {
  my $intervals = shift;

  return [ map { $_->[1] } sort { $a->[0] <=> $b->[0] } map { [ $_->start, $_ ] } @{$intervals} ];
}

sub uniq {
  my $intervals = shift;

  tie my %seen, 'Tie::RefHash';
  
  return [ grep { ! $seen{ $_ }++ } @{$intervals} ];
}

1;

