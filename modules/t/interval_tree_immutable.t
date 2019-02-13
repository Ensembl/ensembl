# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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

throws_ok { Bio::EnsEMBL::Utils::Interval->new() } qr/specify.+?boundaries/, 'Throws with no arguments';
throws_ok { Bio::EnsEMBL::Utils::Interval->new(1) } qr/specify.+?boundaries/, 'Throws with an undefined argument';
throws_ok { Bio::EnsEMBL::Utils::Interval->new(10, 1) } qr/start.+?end/, 'Throws with invalid arguments';
throws_ok { Bio::EnsEMBL::Utils::Interval->new(100, 10) } qr/start.+?end/, 'Throws with invalid arguments';

# degenerate (point) case
my $i = Bio::EnsEMBL::Utils::Interval->new(10, 10);
isa_ok($i, 'Bio::EnsEMBL::Utils::Interval');
ok($i->is_empty, 'empty interval');
ok($i->is_point, 'interval is point');

# a normal interval, start < end
$i = Bio::EnsEMBL::Utils::Interval->new(100, 200);
isa_ok($i, 'Bio::EnsEMBL::Utils::Interval');
is($i->start, 100, 'start position');
is($i->end, 200, 'end position');
ok(!$i->is_empty, 'interval not empty');
ok(!$i->is_point, 'interval\'s not a point');
ok($i->contains(100) && $i->contains(200) && $i->contains(150), 'interval contains points');
ok(!$i->contains(99) && !$i->contains(201), 'interval does not contain points');
# check is_right_of/is_left_of with point/interval
ok(!$i->is_right_of && !$i->is_left_of, 'interval is not left/right of nothing');
ok($i->is_right_of(99), 'interval right of point');
ok(!$i->is_right_of(100) && !$i->is_right_of(150) && !$i->is_right_of(201), 'interval not right of point');
ok($i->is_left_of(201), 'interval left of point');
ok(!$i->is_left_of(99) && !$i->is_left_of(150) && !$i->is_left_of(200), 'interval not left of point');
my $j = Bio::EnsEMBL::Utils::Interval->new(50, 99);
my $k = Bio::EnsEMBL::Utils::Interval->new(50, 150);
my $l = Bio::EnsEMBL::Utils::Interval->new(201, 250);
ok($i->is_right_of($j), 'interval right of another');
ok(!$i->is_right_of($k) && !$i->is_right_of($l), 'interval not right of others');
ok($i->is_left_of($l), 'interval left of another');
ok(!$i->is_left_of($j) && !$i->is_left_of($k), 'interval not left of others');

# check interval data
$j = Bio::EnsEMBL::Utils::Interval->new(100, 200, [100, 200]);
is_deeply($j->data, [100, 200], 'interval data');

# check intersection with other intervals
$k = Bio::EnsEMBL::Utils::Interval->new(50, 150);
ok($i->intersects($k), 'intervals intersect');
$k = Bio::EnsEMBL::Utils::Interval->new(150, 250);
ok($i->intersects($k), 'intervals intersect');
$k = Bio::EnsEMBL::Utils::Interval->new(50, 99);
ok(!$i->intersects($k), 'intervals do not intersect');
$k = Bio::EnsEMBL::Utils::Interval->new(201, 250);
ok(!$i->intersects($k), 'intervals do not intersect');

use_ok 'Bio::EnsEMBL::Utils::Tree::Interval::Immutable::Node';

use_ok 'Bio::EnsEMBL::Utils::Tree::Interval::Immutable';

my $intervals = [ Bio::EnsEMBL::Utils::Interval->new(121626874, 122092717),
		  Bio::EnsEMBL::Utils::Interval->new(121637917, 121658918),
		  Bio::EnsEMBL::Utils::Interval->new(122096077, 124088369) ];
my $tree = Bio::EnsEMBL::Utils::Tree::Interval::Immutable->new($intervals);
isa_ok($tree, 'Bio::EnsEMBL::Utils::Tree::Interval::Immutable');
isa_ok($tree->root, 'Bio::EnsEMBL::Utils::Tree::Interval::Immutable::Node');

# check in-order traversal
my $in_order_results = $tree->in_order_traversal;
is(scalar @{$in_order_results}, 3, 'Number of in-order traversal results');
ok($in_order_results->[0]->start == 121626874 && $in_order_results->[0]->end == 122092717, 'in-order traversal result');
ok($in_order_results->[1]->start == 121637917 && $in_order_results->[1]->end == 121658918, 'in-order traversal result');
ok($in_order_results->[2]->start == 122096077 && $in_order_results->[2]->end == 124088369, 'in-order traversal result');

my ($start, $end, $result);
$result = $tree->query($start, $end);
cmp_deeply($result, [], 'empty query');
$result = $tree->query(undef, $end);
cmp_deeply($result, [], 'empty query');
my $result1 = $tree->query(121779004);
my $result2 = $tree->query(121779004, 121779004);
is(scalar @{$result1}, 1, 'number of query results');
isa_ok($result1->[0], 'Bio::EnsEMBL::Utils::Interval');
ok($result1->[0]->start == 121626874 && $result1->[0]->end == 122092717, 'query result interval boundaries');
cmp_deeply($result1, $result2, 'same result, different query interface');

done_testing();
