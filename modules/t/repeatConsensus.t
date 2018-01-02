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
use Bio::EnsEMBL::Test::TestUtils;

use Bio::EnsEMBL::RepeatConsensus;

use Test::More;
use Test::Warnings;
my $verbose = 0;

my $consensus = 'actg';
my $name   =  'actg(n)';
my $length =  4;
my $class  = 'Simple_repeat';
my $dbID  = 123;

#
# Test constructor
#
my $repeat_consensus = Bio::EnsEMBL::RepeatConsensus->new
  (-REPEAT_CONSENSUS => $consensus,
   -NAME             => $name,
   -LENGTH           => $length,
   -REPEAT_CLASS    => $class,
   -DBID             => 123);

ok ($repeat_consensus && ref($repeat_consensus) && 
    $repeat_consensus->isa('Bio::EnsEMBL::RepeatConsensus'));

ok($repeat_consensus->length() == $length);
ok($repeat_consensus->repeat_consensus() eq $consensus);
ok($repeat_consensus->seq() eq $consensus);
ok($repeat_consensus->name() eq $name);
ok($repeat_consensus->dbID() == $dbID);
ok($repeat_consensus->repeat_class() eq $class);

ok(test_getter_setter($repeat_consensus,'length',10));
ok(test_getter_setter($repeat_consensus,'repeat_class','dummy'));
ok(test_getter_setter($repeat_consensus,'name','dummy'));
ok(test_getter_setter($repeat_consensus,'repeat_consensus','ATGCATGCAT'));
ok(test_getter_setter($repeat_consensus,'dbID',42));

done_testing();
