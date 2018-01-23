# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

use strict;
use warnings;

use Test::More;
use Test::Deep;
use Test::Exception;
use Test::Warnings qw(warning);

use Data::Dumper;

use_ok 'Bio::EnsEMBL::Utils::Tree::Interval::Interval';
use_ok 'Bio::EnsEMBL::Utils::Tree::Interval';

srand(2791);

my $intervals = [ Bio::EnsEMBL::Utils::Tree::Interval::Interval->new(121626874, 122092717),
		  Bio::EnsEMBL::Utils::Tree::Interval::Interval->new(121637917, 121658918),
		  Bio::EnsEMBL::Utils::Tree::Interval::Interval->new(122096077, 124088369) ];

my $tree = Bio::EnsEMBL::Utils::Tree::Interval->new();
isa_ok($tree, 'Bio::EnsEMBL::Utils::Tree::Interval');
ok(!$tree->root, 'Empty tree');
ok(!$tree->search($intervals->[0]), 'Empty tree search');
ok(!$tree->remove($intervals->[0]), 'Empty tree removal');

my $query = make_interval(50, 100, 'data');
$tree->insert($query);
my $search_result = $tree->search($query);
is(scalar @{$search_result}, 1, 'Number of search results');
is($search_result->[0]->data, 'data', 'Correctly insert into empty tree');

$tree = Bio::EnsEMBL::Utils::Tree::Interval->new();
$tree->insert(make_interval(50, 150, 'data1'));
$tree->insert(make_interval(50, 100, 'data2'));
$search_result = $tree->search(make_interval(75, 100));
is(scalar @{$search_result}, 2, 'Number of search results');
is($search_result->[0]->data, 'data1', 'Correct insertion into node with same key');
is($search_result->[1]->data, 'data2', 'Correct insertion into node with same key');
$search_result = $tree->search(make_interval(125, 150));
is(scalar @{$search_result}, 1, 'Number of search results');
is($search_result->[0]->data, 'data1', 'Correct insertion into node with same key');

$tree = Bio::EnsEMBL::Utils::Tree::Interval->new();
$tree->insert(make_interval(50, 150, 'data1'));
$tree->insert(make_interval(25, 100, 'data2'));
# $search_result = $tree->search(make_interval(75, 100));
# is(scalar @{$search_result}, 2, 'Number of search results');
# is($search_result->[0]->data, 'data2', 'Correct insertion into left subtree');
# is($search_result->[1]->data, 'data1', 'Correct insertion into left subtree');

# note Dumper $search_result;

# # check in-order traversal
# my $in_order_results = $tree->in_order_traversal;
# is(scalar @{$in_order_results}, 3, 'Number of in-order traversal results');
# ok($in_order_results->[0]->start == 121626874 && $in_order_results->[0]->end == 122092717, 'in-order traversal result');
# ok($in_order_results->[1]->start == 121637917 && $in_order_results->[1]->end == 121658918, 'in-order traversal result');
# ok($in_order_results->[2]->start == 122096077 && $in_order_results->[2]->end == 124088369, 'in-order traversal result');

# my ($start, $end, $result);
# $result = $tree->query($start, $end);
# cmp_deeply($result, [], 'empty query');
# $result = $tree->query(undef, $end);
# cmp_deeply($result, [], 'empty query');
# my $result1 = $tree->query(121779004);
# my $result2 = $tree->query(121779004, 121779004);
# is(scalar @{$result1}, 1, 'number of query results');
# isa_ok($result1->[0], 'Bio::EnsEMBL::Utils::CenteredIntervalTree::Interval');
# ok($result1->[0]->start == 121626874 && $result1->[0]->end == 122092717, 'query result interval boundaries');
# cmp_deeply($result1, $result2, 'same result, different query interface');

# for my $i (1 .. 100) {
#   my ($start, $end) = (make_random_int(0, 100), make_random_int(0, 100));
  
#   if ($start > $end) {
#     my $tmp = $end;
#     $end = $start;
#     $start = $tmp;
#   }

#   $random_tree->insert(make_interval($start, $end, $data));
# }

sub make_interval {
  my ($start, $end, $data) = @_;

  return Bio::EnsEMBL::Utils::Tree::Interval::Interval->new($start, $end, $data);
}

sub make_random_int {
  my ($min, $max) = @_;

  return $min + rand($max) + 1;
}

done_testing();
