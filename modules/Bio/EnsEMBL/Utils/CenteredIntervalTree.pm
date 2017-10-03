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

sub search {
  my ($self, $start, $end) = @_;
  defined $start or throw 'Undefined point or interval';

  # interval search
  if (defined $end) {
    my $result = [];
    for my $i ($start .. $end + 1) {
      push @{$result}, @{$self->search($i)};
    }

    return sort_by_begin(uniq($result));
  } 

  # point search 
  return $self->_search($self->{top_node}, $start, []);
}

sub _search {
  my ($self, $node, $point, $result) = @_;

  # can be optimized, see https://en.wikipedia.org/wiki/Interval_tree
  for my $i (@{$node->scenter}) {
    push @{$result}, $i if $i->start <= $point and $i->end >= $point;
  }

  if ($point < $node->xcenter and $node->left_child) {
    my $left_results = $self->_search($node->left_child, $point, []);
    push @{$result}, @{$left_results} if $left_results;
  }

  if ($point > $node->xcenter and $node->right_child) {
    my $right_results = $self->_search($node->right_child, $point, []);
    push @{$result}, @{$right_results} if $right_results;
  } 

  return uniq($result);
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

