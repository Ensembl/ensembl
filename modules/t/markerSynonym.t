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
use Bio::EnsEMBL::Map::MarkerSynonym;

use Bio::EnsEMBL::Test::TestUtils;

our $verbose = 0; #set to 1 to turn on debug printouts


######
# 1  #
######

#test constructor
my $source = 'genbank';
my $name = 'DS1234';
my $dbID = 10;


my $ms = Bio::EnsEMBL::Map::MarkerSynonym->new($dbID, $source, $name);

ok($ms && ref $ms && $ms->isa('Bio::EnsEMBL::Map::MarkerSynonym'));

#######
# 2-3 #
#######

#test source

ok($source eq $ms->source);
ok(&test_getter_setter($ms, 'source', 'uniSTS'));

#######
# 4-5 #
#######

#test name
ok($name eq $ms->name);
ok(&test_getter_setter($ms, 'name', 'new_name'));


#######
# 6-7 #
#######

#test dbID
ok($dbID == $ms->dbID);
ok(&test_getter_setter($ms, 'dbID', 123)); 

done_testing();
