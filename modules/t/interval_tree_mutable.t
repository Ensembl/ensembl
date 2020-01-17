# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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

use_ok 'Bio::EnsEMBL::Utils::Interval';
use_ok 'Bio::EnsEMBL::Utils::Tree::Interval::Mutable';

my $intervals = [ Bio::EnsEMBL::Utils::Interval->new(121626874, 122092717),
		  Bio::EnsEMBL::Utils::Interval->new(121637917, 121658918),
		  Bio::EnsEMBL::Utils::Interval->new(122096077, 124088369) ];

my $tree = Bio::EnsEMBL::Utils::Tree::Interval::Mutable->new();
isa_ok($tree, 'Bio::EnsEMBL::Utils::Tree::Interval::Mutable');
ok(!$tree->root, 'Empty tree');
ok(!$tree->search($intervals->[0]), 'Empty tree search');
ok(!$tree->remove($intervals->[0]), 'Empty tree removal');
ok(!$tree->size(), 'Empty tree size');

$tree->insert(make_interval(50, 100, 'data'));
my $search_result = $tree->search(50, 100);
is(scalar @{$search_result}, 1, 'Number of search results');
is($search_result->[0]->data, 'data', 'Correctly insert into empty tree');
is($tree->size(), 1, 'Tree size');

$tree = Bio::EnsEMBL::Utils::Tree::Interval::Mutable->new();
$tree->insert(make_interval(50, 150, 'data1'));
$tree->insert(make_interval(50, 100, 'data2'));
is($tree->size(), 2, 'Tree size');
$search_result = $tree->search(75, 100);
is(scalar @{$search_result}, 2, 'Number of search results');
is($search_result->[0]->data, 'data1', 'Correct insertion into node with same key');
is($search_result->[1]->data, 'data2', 'Correct insertion into node with same key');
$search_result = $tree->search(125, 150);
is(scalar @{$search_result}, 1, 'Number of search results');
is($search_result->[0]->data, 'data1', 'Correct insertion into node with same key');

$tree = Bio::EnsEMBL::Utils::Tree::Interval::Mutable->new();
$tree->insert(make_interval(50, 150, 'data1'));
$tree->insert(make_interval(25, 100, 'data2'));
is($tree->size(), 2, 'Tree size');
$search_result = $tree->search(75, 100);
is(scalar @{$search_result}, 2, 'Number of search results');
is($search_result->[0]->data, 'data2', 'Correct insertion into left subtree');
is($search_result->[1]->data, 'data1', 'Correct insertion into left subtree');

$tree = Bio::EnsEMBL::Utils::Tree::Interval::Mutable->new();
$tree->insert(make_interval(50, 150, 'data1'));
$tree->insert(make_interval(75, 100, 'data2'));
is($tree->size(), 2, 'Tree size');
$search_result = $tree->search(85, 100);
is(scalar @{$search_result}, 2, 'Number of search results');
is($search_result->[0]->data, 'data1', 'Correct insertion into right subtree');
is($search_result->[1]->data, 'data2', 'Correct insertion into right subtree');

$tree = Bio::EnsEMBL::Utils::Tree::Interval::Mutable->new();
$tree->insert(make_interval(50, 150, 'data1'));
$tree->insert(make_interval(75, 200, 'data2'));
is($tree->size(), 2, 'Tree size');
$search_result = $tree->search(50, 100);
is(scalar @{$search_result}, 2, 'Number of search results');
is($search_result->[0]->data, 'data1', 'Correct inclusion');
is($search_result->[1]->data, 'data2', 'Correct inclusion');
$search_result = $tree->search(0, 50);
is(scalar @{$search_result}, 1, 'Number of search results');
is($search_result->[0]->data, 'data1', 'Correct inclusion');
$search_result = $tree->search(200, 300);
is(scalar @{$search_result}, 1, 'Number of search results');
is($search_result->[0]->data, 'data2', 'Correct inclusion');

$tree = Bio::EnsEMBL::Utils::Tree::Interval::Mutable->new();
ok(!$tree->remove(make_interval(50, 100)), 'Delete from empty tree');

$tree = Bio::EnsEMBL::Utils::Tree::Interval::Mutable->new();
$tree->insert(make_interval(75, 150, 'data'));
is($tree->size(), 1, 'Tree size');
ok($tree->remove(make_interval(75, 150)), 'Delete from non-empty tree');
ok(!$tree->size(), 'Tree size');
ok(!$tree->root, 'Root deleted');

$tree = Bio::EnsEMBL::Utils::Tree::Interval::Mutable->new();
$tree->insert(make_interval(50, 120, 'data1'));
$tree->insert(make_interval(75, 100, 'data2'));
$tree->insert(make_interval(75, 200, 'First data to remove'));
$tree->insert(make_interval(75, 150, 'Second data to remove'));
is($tree->size(), 4, 'Tree size');
$search_result = $tree->search(50, 200);
is(scalar @{$search_result}, 4, 'Number of search results');
is($search_result->[0]->data, 'data1', 'Search result');
is($search_result->[1]->data, 'data2', 'Search result');
is($search_result->[2]->data, 'First data to remove', 'Search result');
is($search_result->[3]->data, 'Second data to remove', 'Search result');
ok($tree->remove(make_interval(75, 200)), 'Delete from node with multiple data');
is($tree->size(), 3, 'Tree size');
$search_result = $tree->search(50, 200);
is(scalar @{$search_result}, 3, 'Number of search results');
is($search_result->[0]->data, 'data1', 'Search result');
is($search_result->[1]->data, 'data2', 'Search result');
is($search_result->[2]->data, 'Second data to remove', 'Search result');
ok($tree->remove(make_interval(75, 150)), 'Delete from node with multiple data');
is($tree->size(), 2, 'Tree size');
$search_result = $tree->search(50, 200);
is(scalar @{$search_result}, 2, 'Number of search results');
is($search_result->[0]->data, 'data1', 'Search result');
is($search_result->[1]->data, 'data2', 'Search result');

$tree = Bio::EnsEMBL::Utils::Tree::Interval::Mutable->new();
throws_ok { $tree->insert(make_interval(200, 100, 'spanning_interval')) }
  qr/Cannot insert an interval that spans the origin into a mutable tree/,
  'exception when trying to insert an interval that spans the origin';

$tree = Bio::EnsEMBL::Utils::Tree::Interval::Mutable->new();
map { $tree->insert($_) } @{$intervals};
is($tree->size(), scalar @{$intervals}, 'Tree size');
$search_result = $tree->search(121779004, 121779004);
is(scalar @{$search_result}, 1, 'Number of search results');
ok($search_result->[0]->start == 121626874 && $search_result->[0]->end == 122092717, 'Query result interval boundaries');
$search_result = $tree->search(1e6, 1e7);
is(scalar @{$search_result}, 0, 'Number of search results');

sub make_interval {
  my ($start, $end, $data) = @_;

  return Bio::EnsEMBL::Utils::Interval->new($start, $end, $data);
}

done_testing();
