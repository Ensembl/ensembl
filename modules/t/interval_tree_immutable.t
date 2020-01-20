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

throws_ok { Bio::EnsEMBL::Utils::Interval->new() } qr/specify.+?boundaries/, 'Throws with no arguments';
throws_ok { Bio::EnsEMBL::Utils::Interval->new(1) } qr/specify.+?boundaries/, 'Throws with an undefined argument';

# degenerate (point) case
my $i = Bio::EnsEMBL::Utils::Interval->new(10, 10);
isa_ok($i, 'Bio::EnsEMBL::Utils::Interval');
ok($i->is_empty, 'empty interval');
ok($i->is_point, 'interval is point');

# a normal interval, start < end
$i = Bio::EnsEMBL::Utils::Interval->new(100, 200);
# an interval spanning the origin, start > end
my $i_span = Bio::EnsEMBL::Utils::Interval->new(200,100);

isa_ok($i, 'Bio::EnsEMBL::Utils::Interval');
isa_ok($i_span, 'Bio::EnsEMBL::Utils::Interval');

is($i->spans_origin, 0, 'spans_origin false for non-spanning interval');
is($i_span->spans_origin, 1, 'spans_origin true for spanning interval');

is($i->start, 100, 'start position');
is($i->end, 200, 'end position');
is($i_span->start, 200, 'spanning start position');
is($i_span->end, 100, 'spanning end position');

ok(!$i->is_empty, 'interval not empty');
ok(!$i->is_point, 'interval\'s not a point');
ok(!$i_span->is_empty, 'interval not empty');
ok(!$i_span->is_point, 'interval\'s not a point');

ok($i->contains(100) && $i->contains(200) && $i->contains(150), 'interval contains points');
ok(!$i->contains(99) && !$i->contains(201), 'interval does not contain points');
ok($i_span->contains(100) && $i_span->contains(200) && $i_span->contains(250), 'spanning interval contains points');
ok(!$i_span->contains(101) && !$i_span->contains(199), 'spanning interval does not contain points');

# check is_right_of/is_left_of with point/interval
ok(!$i->is_right_of && !$i->is_left_of, 'interval is not left/right of nothing');
ok(!$i_span->is_right_of && !$i_span->is_left_of, 'spanning interval is not left/right of nothing');

ok($i->is_right_of(99), 'interval right of point');
ok(!$i->is_right_of(100) && !$i->is_right_of(150) && !$i->is_right_of(201), 'interval not right of point');
ok($i->is_left_of(201), 'interval left of point');
ok(!$i->is_left_of(99) && !$i->is_left_of(150) && !$i->is_left_of(200), 'interval not left of point');
throws_ok { $i_span->is_right_of(150) }
  qr/is_right_of not defined for an interval that spans the origin/,
  'exception calling is_right_of with a spanning interval and a point';
throws_ok { $i_span->is_left_of(150) }
  qr/is_left_of not defined for an interval that spans the origin/,
  'exception calling is_left_of with a spanning interval and a point';

my $j = Bio::EnsEMBL::Utils::Interval->new(50, 99);
my $k = Bio::EnsEMBL::Utils::Interval->new(50, 150);
my $l = Bio::EnsEMBL::Utils::Interval->new(201, 250);
my $m = Bio::EnsEMBL::Utils::Interval->new(101, 199);
my $n_span = Bio::EnsEMBL::Utils::Interval->new(201,100);

# non-spanning with non-spanning query
ok($i->is_right_of($j), 'interval right of another');
ok(!$i->is_right_of($k) && !$i->is_right_of($l), 'interval not right of others');
ok($i->is_left_of($l), 'interval left of another');
ok(!$i->is_left_of($j) && !$i->is_left_of($k), 'interval not left of others');

# non-spanning with spanning query
throws_ok { $i->is_right_of($n_span) }
  qr/is_right_of not defined for an interval that spans the origin/,
  'exception calling is_right_of with a spanning interval';
throws_ok { $i->is_left_of($n_span) }
  qr/is_left_of not defined for an interval that spans the origin/,
  'exception calling is_left_of with a spanning interval';

# spanning with non-spanning query
throws_ok { $i_span->is_right_of($m) }
  qr/is_right_of not defined for an interval that spans the origin/,
  'exception calling is_right_of with a spanning interval';
throws_ok { $i_span->is_left_of($m) }
  qr/is_left_of not defined for an interval that spans the origin/,
  'exception calling is_left_of with a spanning interval';

# spanning with spanning query
throws_ok { $i_span->is_right_of($n_span) }
  qr/is_right_of not defined for an interval that spans the origin/,
  'exception calling is_right_of with a spanning interval';
throws_ok { $i_span->is_left_of($n_span) }
  qr/is_left_of not defined for an interval that spans the origin/,
  'exception calling is_left_of with a spanning interval';


# check interval data
$j = Bio::EnsEMBL::Utils::Interval->new(100, 200, [100, 200]);
is_deeply($j->data, [100, 200], 'interval data');

# check intersection with other intervals
$k = Bio::EnsEMBL::Utils::Interval->new(50, 150);
ok($i->intersects($k), 'intervals intersect');
ok($i_span->intersects($k), 'spanning interval and interval intersect');
$k = Bio::EnsEMBL::Utils::Interval->new(150, 250);
ok($i->intersects($k), 'intervals intersect');
$k = Bio::EnsEMBL::Utils::Interval->new(50, 99);
ok(!$i->intersects($k), 'intervals do not intersect');
$k = Bio::EnsEMBL::Utils::Interval->new(201, 250);
ok(!$i->intersects($k), 'intervals do not intersect');
$k = Bio::EnsEMBL::Utils::Interval->new(101,199);
ok(!$i_span->intersects($k), 'spanning interval and interval do not intersect');
ok($i_span->intersects($n_span), 'two spanning intervals intersect');
ok($i->intersects($n_span), 'interval and spanning interval intersect');
my $o_span = Bio::EnsEMBL::Utils::Interval->new(201,99);
ok(!$i->intersects($o_span), 'interval and spanning interval do not intersect');

use_ok 'Bio::EnsEMBL::Utils::Tree::Interval::Immutable::Node';

use_ok 'Bio::EnsEMBL::Utils::Tree::Interval::Immutable';

my $intervals_with_span = [ Bio::EnsEMBL::Utils::Interval->new(20, 30),
                            Bio::EnsEMBL::Utils::Interval->new(30, 20)];

throws_ok { my $impossible_tree = Bio::EnsEMBL::Utils::Tree::Interval::Immutable->new($intervals_with_span) }
  qr/Cannot build a tree containing an interval that spans the origin/,
  'exception when building an interval tree with an interval that spans the origin';

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
