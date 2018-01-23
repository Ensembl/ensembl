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
use Bio::EnsEMBL::Test::TestUtils;

use Bio::EnsEMBL::Map::DBSQL::DitagAdaptor;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db    = $multi->get_DBAdaptor( 'core' );

my $name      = "101A01-2";
my $type      = "ZZ13";
my $tag_count = 2;
my $sequence  = "GAGAACTTGGACCGCAGAGAATACACACAAATCAAACC";

######
# 1  #
######

#test get_DitagAdaptor
my $ditag_adaptor = $db->get_DitagAdaptor;
ok($ditag_adaptor && ref $ditag_adaptor);

#####
# 2 #
#####

#test store
my $new_ditag = Bio::EnsEMBL::Map::Ditag->new (
  -name     => $name, 
  -type     => $type,
	-count    => $tag_count,
  -sequence => $sequence,
);
my @ditags = ( $new_ditag );

#hide the contents of ditag table
$multi->hide('core', 'ditag');

ok($ditag_adaptor->store(\@ditags));

#######
# 3-4 #
#######

#test fetch_all_by_name

my $ditags = $ditag_adaptor->fetch_all_by_name($name);
#if feature was stored it has a dbID now
ok(scalar(@$ditags) && $ditags->[0]->isa('Bio::EnsEMBL::Map::Ditag'));
ok($ditags->[0]->name eq $name);

#######
# 5-6 #
#######

#test fetch_all_by_type

$ditags = $ditag_adaptor->fetch_all_by_type($type);
ok(scalar(@$ditags) && $ditags->[0]->isa('Bio::EnsEMBL::Map::Ditag'));
ok($ditags->[0]->type eq $type);

#######
# 7-8 #
#######

#test fetch_by_name_and_type

my $ditag = $ditag_adaptor->fetch_by_name_and_type($name, $type);
ok($ditag && $ditag->isa('Bio::EnsEMBL::Map::Ditag'));
ok($ditag->type eq $type && $ditag->name eq $name);

######
# 9  #
######

#test fetch_all
$ditags = $ditag_adaptor->fetch_all();
ok(scalar(@$ditags) && $ditags->[0]->isa('Bio::EnsEMBL::Map::Ditag'));

#unhide table after having stored
$multi->restore('core', 'ditag');

######
# 10 #
######

#test fetch_all_with_limit
$ditags = $ditag_adaptor->fetch_all_with_limit($type, 5, 0);
ok((scalar(@$ditags) == 5) && $ditags->[0]->isa('Bio::EnsEMBL::Map::Ditag'));

######
# 11 #
######

#test list_dbIDs

my $dbIDs = $ditag_adaptor->list_dbIDs();
ok(scalar @$dbIDs);

done_testing();
