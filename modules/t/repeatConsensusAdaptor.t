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
use Bio::EnsEMBL::Test::TestUtils;

use Bio::EnsEMBL::Test::MultiTestDB;

my $multi_db = Bio::EnsEMBL::Test::MultiTestDB->new;
my $db = $multi_db->get_DBAdaptor('core');

my $verbose = 0;

# Test Creation

my $rca = $db->get_RepeatConsensusAdaptor();

ok(ref($rca) && $rca->isa('Bio::EnsEMBL::DBSQL::RepeatConsensusAdaptor'));

#
# Test fetch_by_dbID
#

my $rc = $rca->fetch_by_dbID(9);
ok($rc->name() eq 'MIR3');
ok($rc->dbID == 9);
ok($rc->repeat_consensus eq '');
ok($rc->length() == 0);
ok($rc->repeat_class eq 'Type I Transposons/SINE');
ok( $rc->repeat_type eq "ALU" );

#
# Test fetch_by_name
#
$rc = $rca->fetch_by_name('MIR');
ok($rc->name() eq 'MIR');
ok($rc->dbID() == 1);
ok($rc->repeat_consensus eq '');
ok($rc->length() == 0);
ok($rc->repeat_class eq 'Type I Transposons/SINE');

#
# Test fetch_by_name_class
#

$rc = $rca->fetch_by_name_class('MER65A', 'LTRs');
ok($rc->name() eq 'MER65A');
ok($rc->dbID() == 283);
ok($rc->repeat_class eq 'LTRs');
ok($rc->repeat_consensus eq '');
ok($rc->length() == 0);

#
# Test fetch_all_by_class_seq
#
ok(@{$rca->fetch_all_by_class_seq('LTRs', '')} == 38);

#
# Test distinct repeat types retrieval
#

is(@{$rca->fetch_all_repeat_types()}, 3, 'Checking number of repeat types returned');

#
# Test store
#

$multi_db->save('core', 'repeat_consensus', 'meta_coord');

$rc = Bio::EnsEMBL::RepeatConsensus->new
  (-REPEAT_CONSENSUS => 'ACTG',
   -REPEAT_TYPE      => "testtype",
   -NAME             => 'ACTG(n)',
   -LENGTH           => 4,
   -REPEAT_CLASS    => 'Simple_repeat');


$rca->store($rc);

ok($rc->dbID && $rc->adaptor());

$rc = $rca->fetch_by_dbID($rc->dbID);

ok($rc->repeat_consensus eq 'ACTG');
ok($rc->repeat_class eq  'Simple_repeat');
ok($rc->length() == 4);
ok($rc->name eq 'ACTG(n)');
ok($rc->repeat_type() eq "testtype" );

$multi_db->restore('core', 'repeat_consensus');

done_testing();
