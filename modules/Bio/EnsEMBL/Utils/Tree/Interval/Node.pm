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

Bio::EnsEMBL::Utils::Tree::Interval::Node

=head1 SYNOPSIS


=head1 DESCRIPTION

=head1 METHODS

=cut

package Bio::EnsEMBL::Utils::Tree::Interval::Node;

use strict;

use Scalar::Util qw(looks_like_number);
use List::Util qw(max);
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw);

=head1 METHODS 

=head2 new

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($tree, $interval) = @_;
  throw 'Node constructor takes (tree, interval) as arguments'
    unless $tree and $interval;

  my $self = bless({ tree     => $tree,
		     interval => $interval,
		     key      => $interval->start,
		     max      => $interval->end,
		     parent => undef,
		     height => 0,
		     left   => undef,
		     right  => undef }, $class);
  
  return $self;
}

=head2 tree

=cut

sub tree {
  return shift->{tree};
}

=head2 key

=cut

sub key {
  return shift->{key};
}

=head2 parent

=cut

sub parent {
  my $self = shift;
  $self->{parent} = shift if( @_ );
  
  return $self->{parent};
}

=head2 height

=cut

sub height {
  my $self = shift;
  $self->{height} = shift if( @_ );
  
  return $self->{height};
}

=head2 left

=cut

sub left {
  return shift->{left};
}

=head2 right

=cut

sub right {
  return shift->{right};
}

=head2 insert

=cut

sub insert {
  my ($self, $i) = @_;

  if ($i->start < $self->key) {
    # insert into left subtree
    unless (defined $self->left) {
      $self->left = Bio::EnsEMBL::Utils::Tree::Interval::Node->new($self->tree, $i);
      $self->left->parent($self);
    } else {
      $self->left->insert($i);
    }
  } else {
    # insert into right subtree
    unless (defined $self->right) {
      $self->right = Bio::EnsEMBL::Utils::Tree::Interval::Node->new($self->tree, $i);
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

=head1 PRIVATE METHODS

=head2 _update_height

=cut

sub _update_height {
  my $self = shift;

  $self->height(List::Util::max $self->left->height, $self->right->height + 1);
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
  
  if ($left->height - $right->height >= 2) {
    if ($left->left->height >= $left->right->height) {
      # Left-Left case
      $self->_right_rotate;
      $self->_update_max_right_rotate;
    } else {
      # Left-Right case
      $left->_left_rotate;
      $self->_right_rotate;
      $self->_update_max_right_rotate;
    }
  } elsif ($right->height - $left->height >= 2) {
    if ($right->right->height >= $right->left->height) {
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

=head2 _left_rotate

=cut

sub _left_rotate {
  my $self = shift;

  my $right = $self->right;

  $right->parent($self->parent);
  unless (defined $right->parent) {
    $self->tree->root($right);
  } else {
    if ($right->parent->left == $self) {
      $right->parent->left = $right
    } elsif ($right->parent->right == $self) {
      $right->parent->right = $right;
    }
  }
  
  $self->right = $right->left;
  $self->right->parent($self) if $self->right;
  
  $right->left = $self;
  $self->parent($right);

  $self->_update_height;
  $right->_update_height;
  
}

=head2 _update_max_left_rotate

=cut

sub _update_max_left_rotate {
  my $self = shift;
}

=head2 _right_rotate

=cut

sub _right_rotate {
  my $self = shift;

  my $left = $self->left;

  $left->parent($self->parent);
  unless (defined $left->parent) {
    $self->tree->root($left);
  } else {
    if ($left->parent->left == $self) {
      $left->parent->left = $left;
    } elsif ($left->parent->right == $self) {
      $left->parent->right = $left;
    }
  }
  
  $self->left = $left->right;
  $self->left->parent($self) if $self->left;
  
  $left->right = $self;
  $self->parent($left);

  $self->_update_height;
  $left->_update_height;
}

=head2 _update_max_right_rotate

=cut

sub _update_max_right_rotate {
  my $self = shift;
}

1;

