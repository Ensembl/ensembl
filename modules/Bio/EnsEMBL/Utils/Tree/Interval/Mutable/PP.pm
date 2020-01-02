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

Bio::EnsEMBL::Utils::Tree::Interval::Mutable::PP

=head1 DESCRIPTION

Pure Perl fall back implementation of a mutable interval tree, uses
augmented AVL binary balanced trees.
 
=head1 METHODS

=cut

package Bio::EnsEMBL::Utils::Tree::Interval::Mutable::PP;

use strict;
use Carp;

use Bio::EnsEMBL::Utils::Tree::Interval::Mutable::Node;
use Bio::EnsEMBL::Utils::Interval;
use Bio::EnsEMBL::Utils::Exception qw(throw);

=head2 new

  Arg []      : none
  Example     : my $tree = Bio::EnsEMBL::Utils::Tree::Mutable();
  Description : Constructor. Creates a new mutable tree instance
  Returntype  : Bio::EnsEMBL::Utils::Tree::Interval::Mutable::PP
  Exceptions  : none
  Caller      : general

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  # mmhh, should probably return just the hash ref
  return bless({ _root => undef, _size => 0 }, $class);
}

=head2 root

  Arg []      : none
  Example     : my $root = $tree->root;
  Description : Returns the tree top node
  Returntype  : Bio::EnsEMBL::Utils::Tree::Interval::Mutable::Node
  Exceptions  : none
  Caller      : general

=cut

sub root {
  my $self = shift;
  $self->{_root} = shift if( @_ );
  
  return $self->{_root};
}

=head2 size

  Arg []      : none
  Example     : print "Tree size is ", $tree->size, "\n";
  Description : Return the size of the tree, i.e. the number of nodes
  Returntype  : scalar
  Exceptions  : none
  Caller      : general

=cut

sub size {
  return shift->{_size};
}

=head2 insert

  Arg [1]     : Bio::EnsEMBL::Utils::Interval
                Interval to insert
  Example     : $tree->insert(Bio::EnsEMBL::Utils::Interval->new(10, 20, 'data'));
  Description : Insert an interval in the tree
  Returntype  : scalar (1), upon success
  Exceptions  : thrown if Interval spans origin (has start > end)
  Caller      : general

=cut

sub insert {
  my ($self, $i) = @_;
  
  if ($i->spans_origin) {
    throw "Cannot insert an interval that spans the origin into a mutable tree";
  }
  # base case: empty tree, assign new node to root
  unless (defined $self->root) {
    $self->root(Bio::EnsEMBL::Utils::Tree::Interval::Mutable::Node->new($self, $i));
    $self->{_size}++;
    
    return 1;
  }

  # check if node exists with the same key
  my $node = $self->root->search_by_key($i->start);
  if ($node) {
    # check the node's intervals if there's already one
    # which is the same as the interval we're trying to insert
    map { return 0 if $_->start == $i->start and $_->end == $i->end } @{$node->intervals};

    # add the interval to the node
    $node->add_interval($i);

    # update max of node and its ancestors if necessary
    if ($node->max < $i->end) {
      $node->max($i->end);
      $node->_update_parents_max if $node->parent;
    }

    $self->{_size}++;
    
    return 1;
  } else {
    # node with the interval's key doesn't exist
    # insert from root node
    $self->root->insert($i);
    $self->{_size}++;
    
    return 1;
  }

  croak "Shouldn't be here";
}

=head2 search

  Arg [1]     : scalar, $start
                the starting point of the interval to search
  Arg [2]     : scalar, $end
                the end point of the interval to search
  Example     : my $result = $tree->search(85, 100);
  Description : Search the intervals in the tree overlapping the query
  Returntype  : Arrayref of Bio::EnsEMBL::Utils::Interval instances
  Exceptions  : none
  Caller      : general

=cut

sub search {
  my ($self, $start, $end) = @_;
  
  return unless $self->root; # empty tree
  return $self->root->search(Bio::EnsEMBL::Utils::Interval->new($start, $end));
  
}

=head2 remove

  Arg [1]     : Bio::EnsEMBL::Utils::Interval
                The interval to remove from the tree
  Example     : $tree->remove($i);
  Description : Remove an interval from the tree
  Returntype  : scalar, 1 if removal is successful, 0 otherwise
  Exceptions  : none
  Caller      : general

=cut

sub remove {
  my ($self, $i) = @_;

  return 0 unless $self->root; # empty tree

  my $node = $self->root->search_by_key($i->start);
  return 0 unless $node;
  
  my $node_intervals = $node->intervals;
  if (scalar @{$node_intervals} > 1) {
    # node with this key has more than this interval. Find it and remove
    my $removed = 0;
    foreach my $j (0 .. $#{$node_intervals}) {
      if ($node_intervals->[$j]->start == $i->start and $node_intervals->[$j]->end == $i->end) {
	splice @{$node_intervals}, $j, 1;
	$removed = 1;
	last;
      }
    }

    if ($removed) {
      $removed = 0;

      # update max of node and its ancestors, if necessary
      if ($i->end == $node->max) {
	my $node_highest_end = $node->_highest_end;
	if ($node->left and $node->right) {
	  $node->{max} = max $node->left->max, $node->rigth->max, $node_highest_end;
	} elsif ($node->left and not $node->right) {
	  $node->{max} = max $node->left->max, $node_highest_end;
	} elsif ($node->right and not $node->left) {
	  $node->{max} = max $node->right->max, $node_highest_end;
	} else {
	  $node->{max} = $node_highest_end;
	}

	$node->parent->_update_parents_max if $node->parent;
      }

      $self->{_size}--;
      
      return 1;
    } else {
      return 0;
    }
  } elsif (scalar @{$node_intervals}) {
    # node with this key has only this interval
    # check if the remaining node's interval is the one we want to remove
    if ($node_intervals->[0]->start == $i->start and $node_intervals->[0]->end == $i->end) {
      # remove the whole node from the tree
      if ($node->key == $self->root->key) {
	# we're removing the root node
	# create dummy node temporarily taking root's parent node
	my $root_parent =
	  Bio::EnsEMBL::Utils::Tree::Interval::Mutable::Node->new($self,
								  Bio::EnsEMBL::Utils::Interval->new($i->start, $i->end));
	$root_parent->left($self->root);
	$self->root->parent($root_parent);
	
	my $removed_node = $self->root->remove($node);
	$self->root($root_parent->left);

	$self->root->parent(undef) if $self->root;
	if ($removed_node) {
	  $removed_node = undef;
	  $self->{_size}--;
	  
	  return 1;
	} else {
	  return 0;
	}
      } else {
	my $removed_node = $self->root->remove($node);
	if ($removed_node) {
	  $removed_node = undef;
	  $self->{_size}--;
	  
	  return 1;
	} else {
	  return 0;
	}
      }
    } else {
      # the remaining record is not the one we want to remove
      return 0;
    }
    
  } else {
    # no records at all in this node, shouldn't happen
    croak "This should not happen";
  }
}

sub _in_order_traversal_delete {
  my ($self, $node) = @_;

  return unless $node;

  $self->_in_order_traversal_delete($node->left);
  $node->left(undef);
  $self->_in_order_traversal_delete($node->right);
  $node->right(undef);

  $node->parent(undef);

}

sub DESTROY {
  my $self = shift;

  $self->_in_order_traversal_delete($self->root);
}

1;

