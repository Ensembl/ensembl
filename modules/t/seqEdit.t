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


use Test::More;
use Test::Warnings;
use strict;
use warnings;

use Bio::EnsEMBL::SeqEdit;
use Bio::EnsEMBL::Attribute;
use Bio::EnsEMBL::Test::TestUtils;

my $code  = 'rna edit';
my $desc  = 'Post transcriptional RNA edit';
my $value = '2 3 ACTG';
my $name  = 'RNA Edit';

my $a = Bio::EnsEMBL::Attribute->new
  (-CODE  => $code,
   -DESCRIPTION  => $desc,
   -VALUE => $value,
   -NAME  => $name);

my $se = Bio::EnsEMBL::SeqEdit->new(-ATTRIB => $a);

ok(ref($se) && $se->isa('Bio::EnsEMBL::SeqEdit'));

ok($se->name() eq $name);
ok($se->description() eq $desc);
ok($se->code() eq $code);

ok($se->alt_seq() eq 'ACTG');
ok($se->start()   == 2);
ok($se->end()     == 3);

ok(test_getter_setter($se, 'start', 10));
ok(test_getter_setter($se, 'end', 12));
ok(test_getter_setter($se, 'alt_seq', 'GGAAA'));
ok(test_getter_setter($se, 'name', 'test name'));
ok(test_getter_setter($se, 'description', 'test desc'));
ok(test_getter_setter($se, 'code', 'test cpde'));

ok($se->length_diff == 2);

my $seq = 'CCCC';
$se->apply_edit(\$seq);
ok($seq eq 'CACTGC');

# test insert before first base
$seq = 'ACTG';
$se->alt_seq('CC');
$se->start(1);
$se->end(0);
ok($se->apply_edit(\$seq));
ok($seq eq 'CCACTG');

ok($se->length_diff() == 2);

# test insert after last base
$seq = 'ACTG';
$se->alt_seq('CC');
$se->start(5);
$se->end(4);
$se->apply_edit(\$seq);
ok($seq eq 'ACTGCC');

ok($se->length_diff() == 2);

# test deletion of entire sequence
$seq = 'ACTG';
$se->alt_seq('');
$se->start(1);
$se->end(4);
$se->apply_edit(\$seq);
ok($seq eq '');

ok($se->length_diff() == -4);

# test replacement of some sequence
$seq = 'ACTG';
$se->alt_seq('TC');
$se->start(2);
$se->end(3);
$se->apply_edit(\$seq);
ok($seq eq 'ATCG');

# test conversion to attribute
$a = $se->get_Attribute();
ok($a->name() eq $se->name());
ok($a->description() eq $se->description());
ok($a->code() eq $se->code());
ok($a->value eq '2 3 TC');

done_testing();
