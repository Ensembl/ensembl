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

Bio::EnsEMBL::Utils::Tree::Interval::Mutable::Node

=head1 DESCRIPTION

Represents a node in the mutable interval tree pure perl implementation.

=head1 METHODS

=cut

package Bio::EnsEMBL::Utils::Tree::Interval::Mutable::Node;

use strict;

use Scalar::Util qw(looks_like_number weaken);
use List::Util qw(max);
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw);

=head1 METHODS 

=head2 new

  Arg [1]     : Bio::EnsEMBL::Utils::Tree::Interval::Mutable::PP
                The tree to which the node belongs
  Arg [2]     : Bio::EnsEMBL::Utils::Interval
  Description : Constructor. Creates a new mutable tree instance node
                associated with the given interval
  Returntype  : Bio::EnsEMBL::Utils::Tree::Interval::Mutable::Node
  Exceptions  : none
  Caller      : general

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($tree, $interval) = @_;
  throw 'Node constructor takes (tree, interval) as arguments'
    unless $tree and $interval;

  my $self = bless({ tree      => $tree,
		     intervals => [ $interval ], # the array of all records with the same key
		     key       => $interval->start,
		     max       => $interval->end,
		     parent    => undef,
		     height    => 0,
		     left      => undef,
		     right     => undef }, $class);
  
  return $self;
}

=head2 tree

  Arg []      : none
  Example     : my $tree = $node->root;
  Description : Returns the tree to which the node belongs
  Returntype  : Bio::EnsEMBL::Utils::Tree::Interval::Mutable::PP
  Exceptions  : none
  Caller      : general

=cut

sub tree {
  return shift->{tree};
}

=head2 key

  Arg []      : none
  Example     : my $key = $node->key;
  Description : Returns the key associated with the node
  Returntype  : scalar
  Exceptions  : none
  Caller      : general

=cut

sub key {
  my $self = shift;
  $self->{key} = shift if( @_ );
  
  return $self->{key};
}

=head2 intervals

  Arg []      : none
  Example     : my $intervals = $node->intervals;
  Description : Returns the intervals associated with the node
  Returntype  : Arrayref of Bio::EnsEMBL::Utils::Interval
  Exceptions  : none
  Caller      : general

=cut

sub intervals {
  return shift->{intervals};
}

=head add_interval 

  Arg []      : none
  Description : Add an interval to the node's set of intervals
  Returntype  : none
  Exceptions  : none
  Caller      : general

=cut

sub add_interval {
  push @{shift->{intervals}}, shift;
}

=head2 parent

  Arg []      : none
  Description : Return the parent of the node in the tree
  Returntype  : none
  Exceptions  : none
  Caller      : general

=cut

sub parent {
  my $self = shift;
  if (@_) {
    $self->{parent} = shift;
    weaken($self->{parent});
  } 
  
  return $self->{parent};
}

=head2 height

  Arg []      : none
  Description : Return the height of the node
  Returntype  : scalar, positive or 0
  Exceptions  : none
  Caller      : general

=cut

sub height {
  my $self = shift;
  $self->{height} = shift if( @_ );
  
  return $self->{height};
}

=head2 left

  Arg []      : none
  Description : Return the node's left child
  Returntype  : Bio::EnsEMBL::Utils::Tree::Interval::Mutable::Node
  Exceptions  : none
  Caller      : general

=cut

sub left {
  my $self = shift;
  $self->{left} = shift if( @_ );
  
  return $self->{left};
}

=head2 right

  Arg []      : none
  Description : Return the node's right child
  Returntype  : Bio::EnsEMBL::Utils::Tree::Interval::Mutable::Node
  Exceptions  : none
  Caller      : general

=cut

sub right {
  my $self = shift;
  $self->{right} = shift if( @_ );
  
  return $self->{right};
}

=head2 search

  Arg [1]     : Bio::EnsEMBL::Utils::Interval
                The interval to search for overlaps in the tree
  Description : Search the node and its successors for overlapping
                intervals with the query
  Returntype  : Arrayref of Bio::EnsEMBL::Utils::Interval
  Exceptions  : none
  Caller      : general

=cut

sub search {
  my ($self, $i) = @_;

  # if interval is to the right of the rightmost point of any interval in this node and
  # all its children, there won't be any matches
  return [] if $i->start > $self->{max};

  my $results = [];
  
  # search left subtree
  if ($self->left and $self->left->{max} >= $i->start) {
    push @{$results}, @{$self->left->search($i)};
  }

  # search this node
  push @{$results}, @{$self->_overlapping_intervals($i)};

  # if interval is to the left of the start of this interval, then
  # it can't be in any child to the right
  return $results if $i->end < $self->key;

  # search right subtree
  push @{$results}, @{$self->right->search($i)} if $self->right;

  return $results;
}

=head2 search_by_key

  Arg [1]     : scalar, $key
                The key to search the node for
  Description : Searches for a node by a 'key' value
  Returntype  : Bio::EnsEMBL::Utils::Tree::Interval::Mutable::Node
  Exceptions  : none
  Caller      : general

=cut

sub search_by_key {
  my ($self, $key) = @_;

  if ($self->key == $key) {
    return $self;
  } elsif ($key < $self->key) {
    return $self->left->search_by_key($key) if $self->left;
  } else {
    return $self->right->search_by_key($key) if $self->right;
  }

}

=head2 insert

  Arg [1]     : Bio::EnsEMBL::Utils::Interval
                The interval to insert
  Description : Insert an interval in the node or in its successors
  Returntype  : none
  Exceptions  : none
  Caller      : general

=cut

sub insert {
  my ($self, $i) = @_;

  if ($i->start < $self->key) {
    # insert into left subtree
    unless (defined $self->left) {
      $self->left(Bio::EnsEMBL::Utils::Tree::Interval::Mutable::Node->new($self->tree, $i));
      $self->left->parent($self);
    } else {
      $self->left->insert($i);
    }
  } else {
    # insert into right subtree
    unless (defined $self->right) {
      $self->right(Bio::EnsEMBL::Utils::Tree::Interval::Mutable::Node->new($self->tree, $i));
      $self->right->parent($self);
    } else {
      $self->right->insert($i);
    }
  }

  # update max value if needed
  $self->{max} = $i->end if $self->{max} < $i->end;

  # update node's height
  $self->_update_height;

  # rebalance to ensure O(logn) time operations
  $self->_rebalance;
}

=head2 remove

  Arg [1]     : Bio::EnsEMBL::Utils::Tree::Interval::Mutable::Node
                The node to remove from the tree
  Description : Remove a node from the tree
  Returntype  : none
  Exceptions  : none
  Caller      : general

=cut

sub remove {
  my ($self, $node) = @_;

  return unless $node;

  my $parent = $self->parent;

  if ($node->key < $self->key) {
    # node to be removed is in left subtree
    return $self->left->remove($node) if $self->left;
    return;
  } elsif ($node->key > $self->key) {
    # node to be removed is in right subtree
    return $self->right->remove($node) if $self->right;
    return;
  } else {
    if ($self->left and $self->right) {
      # node has two children
      my $lowest = $self->right->_lowest;
      $self->key($lowest->key);
      $self->{intervals} = $lowest->{intervals};
      return $self->right->remove($self);
    } elsif ($parent->left == $self) {
      # one child or no child case on left side
      if ($self->right) {
	$parent->left($self->right);
	$self->right->parent($parent);
      } else {
	$parent->left($self->left);
	$self->left->parent($parent) if $self->left
      }
      $parent->_update_parents_max;
      $parent->_update_height;
      $parent->_rebalance;

      return $self;
      
    } elsif ($parent->right == $self) {
      # one child or no child case on right side
      if ($self->right) {
	$parent->right($self->right);
	$self->right->parent($parent);
      } else {
	$parent->right($self->left);
	$self->left->parent($parent) if $self->left;
      }
      $parent->_update_parents_max;
      $parent->_update_height;
      $parent->_rebalance;

      return $self;
    }
  }
}

=head1 PRIVATE METHODS

=head2 _height

Not a method since code could invoke method on undefined instances, e.g. _rebalance

=cut

sub _height {
  my $node = shift;

  return -1 unless $node;
  return $node->height;
}

=head2 _lowest 

Returns the 'smallest' node in the tree

=cut

sub _lowest {
  my $self = shift;

  return $self unless $self->left;
  return $self->left->_lowest;
}

sub _highest_end {
  my $self = shift;

  my $high = $self->{intervals}[0]->end;
  map { $high = $_->end if $high < $_->end } @{$self->{intervals}};

  return $high;
}

sub _update_height {
  my $self = shift;

  $self->height(List::Util::max $self->left?$self->left->height:0, $self->right?$self->right->height:0 + 1);
}

sub _update_parents_max {
  my $self = shift;
  # updates the max value of all the parents after inserting into already existing node, as well as
  # removing the node completely or removing the record of an already existing node. Starts with
  # the parent of an affected node and bubbles up to root
  my $high = $self->_highest_end;
  if ($self->left and $self->right) {
    $self->{max} = max $self->left->{max}, $self->right-{max}, $high;
  } elsif ($self->left and !$self->right) {
    $self->{max} = max $self->left->{max}, $high;
  } elsif (!$self->left and $self->right) {
    $self->{max} = max $self->right->{max}, $high;
  } else {
    $self->{max} = $high;
  }

  $self->parent->_update_parents_max if $self->parent;
}

=head2 _rebalance

Rebalances the tree if the height value between two nodes of the same parent is greater than two. 
There are 4 cases that can happen:

Left-Left case:
         z                                      y
        / \                                   /   \
       y   T4      Right Rotate (z)          x     z
      / \          - - - - - - - - ->       / \   / \
     x   T3                                T1 T2 T3 T4
    / \
  T1   T2

Left-Right case:
       z                               z                           x
      / \                             / \                        /   \
     y   T4  Left Rotate (y)         x  T4  Right Rotate(z)     y     z
    / \      - - - - - - - - ->     / \      - - - - - - - ->  / \   / \
  T1   x                           y  T3                      T1 T2 T3 T4
      / \                         / \
    T2   T3                      T1 T2 

Right-Right case:
    z                               y
   / \                            /   \
  T1  y     Left Rotate(z)       z     x
     / \   - - - - - - - ->     / \   / \
    T2  x                      T1 T2 T3 T4
       / \
      T3 T4

Right-Left case:
     z                            z                            x
    / \                          / \                         /   \
   T1  y   Right Rotate (y)     T1  x      Left Rotate(z)   z     y
      / \  - - - - - - - - ->      / \   - - - - - - - ->  / \   / \
     x  T4                        T2  y                   T1 T2 T3 T4
    / \                              / \
  T2   T3                           T3 T4

=cut 

sub _rebalance {
  my $self = shift;

  my ($left, $right) = ($self->left, $self->right);
  
  if (_height($left) - _height($right) >= 2) {
    if (_height($left->left) >= _height($left->right)) {
      # Left-Left case
      $self->_right_rotate;
      $self->_update_max_right_rotate;
    } else {
      # Left-Right case
      $left->_left_rotate;
      $self->_right_rotate;
      $self->_update_max_right_rotate;
    }
  } elsif (_height($right) - _height($left) >= 2) {
    if (_height($right->right) >= _height($right->left)) {
      # Right-Right case
      $self->_left_rotate;
      $self->_update_max_left_rotate;
    } else {
      # Right-Left case
      $right->_right_rotate;
      $self->_left_rotate;
      $self->_update_max_left_rotate;
    }
  }
}

sub _left_rotate {
  my $self = shift;

  my $right = $self->right;

  $right->parent($self->parent);
  unless (defined $right->parent) {
    $self->tree->root($right);
  } else {
    if ($right->parent->left == $self) {
      $right->parent->left($right);
    } elsif ($right->parent->right == $self) {
      $right->parent->right($right);
    }
  }
  
  $self->right($right->left);
  $self->right->parent($self) if $self->right;
  
  $right->left($self);
  $self->parent($right);

  $self->_update_height;
  $right->_update_height;
  
}

sub _update_max_left_rotate { # handles Right-Right case and Right-Left case in rebalancing AVL tree
  my $self = shift;

  # update max of left sibling (x in first case, y in second)
  my $parent = $self->parent;
  my $parent_right = $parent->right;
  if ($parent_right) {
    my $parent_right_high = $parent_right->_highest_end;
    if (!$parent_right->left and $parent_right->right) {
      $parent_right->{max} = max $parent_right_high, $parent_right->right->{max};
    } elsif ($parent_right->left and !$parent_right->right) {
      $parent_right->{max} = max $parent_right_high, $parent_right->left->{max};
    } elsif (!$parent_right->left and !$parent_right->right) {
      $parent_right->{max} = $parent_right_high;
    } else {
      $parent_right->{max} = max $parent_right_high, $parent_right->left->{max}, $parent_right->right->{max};
    }
  }
  
  # update max of itself (z)
  my $high = $self->_highest_end;
  if (!$self->left and $self->right) {
    $self->{max} = max $high, $self->right->{max};
  } elsif ($self->left and !$self->right) {
    $self->{max} = max $high, $self->left->{max};
  } elsif (!$self->left and !$self->right) {
    $self->{max} = $high;
  } else {
    $self->{max} = max $high, $self->left->{max}, $self->right->{max};
  }
  
  # update max of parent (y in first case, x in second)
  $parent->{max} = max $parent->left?$parent->left->{max}:0, $parent->right?$parent->right->{max}:0, $parent->_highest_end
    if $parent;

}

sub _right_rotate {
  my $self = shift;

  my $left = $self->left;

  $left->parent($self->parent);
  unless (defined $left->parent) {
    $self->tree->root($left);
  } else {
    if ($left->parent->left == $self) {
      $left->parent->left($left);
    } elsif ($left->parent->right == $self) {
      $left->parent->right($left);
    }
  }
  
  $self->left($left->right);
  $self->left->parent($self) if $self->left;
  
  $left->right($self);
  $self->parent($left);

  $self->_update_height;
  $left->_update_height;
}

sub _update_max_right_rotate { # handles Left-Left case and Left-Right case after rebalancing AVL tree
  my $self = shift;

  # update max of left sibling (x in first case, y in second)
  my $parent = $self->parent;
  my $parent_left = $parent->left;
  if ($parent_left) {
    my $parent_left_high = $parent_left->_highest_end;
    if (!$parent_left->left and $parent_left->right) {
      $parent_left->{max} = max $parent_left_high, $parent_left->right->{max};
    } elsif ($parent_left->left and !$parent_left->right) {
      $parent_left->{max} = max $parent_left_high, $parent_left->left->{max};
    } elsif (!$parent_left->left and !$parent_left->right) {
      $parent_left->{max} = $parent_left_high;
    } else {
      $parent_left->{max} = max $parent_left_high, $parent_left->left->{max}, $parent_left->right->{max};
    }
  }
  
  # update max of itself (z)
  my $high = $self->_highest_end;
  if (!$self->left and $self->right) {
    $self->{max} = max $high, $self->right->{max};
  } elsif ($self->left and !$self->right) {
    $self->{max} = max $high, $self->left->{max};
  } elsif (!$self->left and !$self->right) {
    $self->{max} = $high;
  } else {
    $self->{max} = max $high, $self->left->{max}, $self->right->{max};
  }
  
  # update max of parent (y in first case, x in second)
  $parent->{max} = max $parent->left?$parent->left->{max}:0, $parent->right?$parent->right->{max}:0, $parent->_highest_end
    if $parent;
}

sub _overlapping_intervals {
  my ($self, $i) = @_;

  my $results = [];
  if ($self->key <= $i->end and $i->start <= $self->_highest_end) {
    # node's interval overlap: check individual intervals
    map { push @{$results}, $_ if $i->start <= $_->end } @{$self->{intervals}}
  }
  
  return $results;
}

1;

