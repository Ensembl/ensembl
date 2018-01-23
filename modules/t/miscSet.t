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

use Test::More;
use Test::Warnings;

use Bio::EnsEMBL::MiscSet;
use Bio::EnsEMBL::Test::TestUtils;

our $verbose = 0; #set to 1 to turn on debug printouts



my $misc_set_id = 3;
my $code = '1mbcloneset';
my $name = '1mb Clone Set';
my $desc = 'This is a 1MB cloneset';
my $max_len = 1e7;


#
# Test constructor
#
my $ms = Bio::EnsEMBL::MiscSet->new('-DBID' => $misc_set_id,
                                    '-CODE' => $code,
                                    '-NAME' => $name,
                                    '-DESCRIPTION' => $desc,
                                    '-LONGEST_FEATURE' => $max_len);

ok($ms->dbID == $misc_set_id);
ok($ms->code eq $code);
ok($ms->description eq $desc);
ok($ms->longest_feature == $max_len);
ok($ms->name eq $name);

ok(test_getter_setter($ms, 'dbID', 12));
ok(test_getter_setter($ms, 'code', 'testcode'));
ok(test_getter_setter($ms, 'description', 'new description'));
ok(test_getter_setter($ms, 'longest_feature', 1e8));
ok(test_getter_setter($ms, 'name', $name));


done_testing();
