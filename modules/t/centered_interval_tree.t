# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2017] EMBL-European Bioinformatics Institute
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
use Test::Exception;
use Test::Warnings qw(warning);

use Data::Dumper;

use_ok 'Bio::EnsEMBL::Utils::CenteredIntervalTree::Interval';

throws_ok { Bio::EnsEMBL::Utils::CenteredIntervalTree::Interval->new() } qr/specify.+?boundaries/, 'Throws with no arguments';
throws_ok { Bio::EnsEMBL::Utils::CenteredIntervalTree::Interval->new(1) } qr/specify.+?boundaries/, 'Throws with an undefined argument';
throws_ok { Bio::EnsEMBL::Utils::CenteredIntervalTree::Interval->new(10, 1) } qr/start.+?end/, 'Throws with invalid arguments';
throws_ok { Bio::EnsEMBL::Utils::CenteredIntervalTree::Interval->new(10, 10) } qr/start.+?end/, 'Throws with invalid arguments';

my $i = Bio::EnsEMBL::Utils::CenteredIntervalTree::Interval->new(1, 10);
isa_ok($i, 'Bio::EnsEMBL::Utils::CenteredIntervalTree::Interval');
is($i->start, 1, 'start position');
is($i->end, 10, 'end position');

my $j = Bio::EnsEMBL::Utils::CenteredIntervalTree::Interval->new(100, 200, [100, 200]);
is_deeply($j->data, [100, 200], 'interval data');

use_ok 'Bio::EnsEMBL::Utils::CenteredIntervalTree::Node';

use_ok 'Bio::EnsEMBL::Utils::CenteredIntervalTree';

my $intervals = [ Bio::EnsEMBL::Utils::CenteredIntervalTree::Interval->new(121626874, 122092717),
		  Bio::EnsEMBL::Utils::CenteredIntervalTree::Interval->new(121637917, 121658918),
		  Bio::EnsEMBL::Utils::CenteredIntervalTree::Interval->new(122096077, 124088369) ];
my $tree = Bio::EnsEMBL::Utils::CenteredIntervalTree->new($intervals);
isa_ok($tree, 'Bio::EnsEMBL::Utils::CenteredIntervalTree');

# print Dumper $tree->search(121779004);

# $it->insert(121626874, 122092717, 'a');
# $it->insert(121637917, 121658918, 'b');
# $it->insert(122096077, 124088369, 'c');

done_testing();
