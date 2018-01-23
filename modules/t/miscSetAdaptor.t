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

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::MiscSet;
use Bio::EnsEMBL::Test::TestUtils;
our $verbose = 0; #set to 1 to turn on debug printouts


my $multi_db = Bio::EnsEMBL::Test::MultiTestDB->new;
my $db = $multi_db->get_DBAdaptor('core');


#
# Test constructor
#

my $msa = $db->get_MiscSetAdaptor();


ok($msa && ref($msa) && $msa->isa('Bio::EnsEMBL::DBSQL::MiscSetAdaptor'));


#
# Test fetch_all
#

ok(@{$msa->fetch_all()} == 4);

#
# Test fetch_by_dbID
#

my $ms = $msa->fetch_by_dbID(1);

ok($ms->dbID(1));
ok($ms->code() eq 'ntctgs');
ok($ms->name() eq 'NT contigs');
ok($ms->description eq '');
ok($ms->longest_feature == 7e7);

#
# Test fetch_by_code
#

$ms = $msa->fetch_by_code('cloneset');
ok($ms->dbID() == 4);
ok($ms->code() eq 'cloneset');
ok($ms->name() eq '1Mb cloneset');
ok($ms->description eq '');
ok($ms->longest_feature == 6e6);


#
# Test store method
#

$multi_db->hide('core', 'misc_set');

my $misc_set = Bio::EnsEMBL::MiscSet->new
  (-CODE => 'code',
   -NAME => 'name',
   -DESCRIPTION => 'description',
   -LONGEST_FEATURE => 1000);


$msa->store($misc_set);

ok($misc_set->dbID());
ok($misc_set->adaptor() == $msa);

my $count = $db->dbc->db_handle->selectall_arrayref
  ("SELECT COUNT(*) FROM misc_set WHERE code = 'code'")->[0]->[0];

ok($count == 1);

# flush and reload cache...
$msa->fetch_all();

$misc_set = $msa->fetch_by_code('code');
ok($misc_set->code() eq 'code');
ok($misc_set->name() eq 'name');
ok($misc_set->description eq 'description');
ok($misc_set->longest_feature == 1000);


# try to store a misc_set with the same code
$misc_set = Bio::EnsEMBL::MiscSet->new
  (-code => 'code',
   -name => 'name',
   -description => 'description',
   -longest_feature => 1000);

$msa->store($misc_set);

ok($misc_set->dbID && $misc_set->adaptor);

$count = $db->dbc->db_handle->selectall_arrayref
  ("SELECT COUNT(*) FROM misc_set WHERE code = 'code'")->[0]->[0];
ok($count == 1);

$multi_db->restore('core', 'misc_set');

done_testing();
