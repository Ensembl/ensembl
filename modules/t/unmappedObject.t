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
use Bio::EnsEMBL::UnmappedObject;
use Bio::EnsEMBL::DBSQL::UnmappedObjectAdaptor;
our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
ok( $multi );

my $db = $multi->get_DBAdaptor( "core" );

my $uma = $db->get_UnmappedObjectAdaptor();

ok(ref($uma));

my $analysis_ad = $db->get_AnalysisAdaptor();

ok($analysis_ad);


#
# Test store
#

$multi->save('core', 'unmapped_object', 'unmapped_reason', 'meta_coord');

my $analysis = $analysis_ad->fetch_by_logic_name('Unigene');


my $test_type = 'xref';
my $test_ex_db_id = 4100;
my $test_identifier = 'X90';
my $test_query_score = 1;
my $test_target_score = 1000;
my $test_summary = 'TEST failed match';
my $test_desc = 'TEST failed match due to being below threshold of 90%';
my $test_ensembl_id = "21734";
my $test_ensembl_object_type = "Translation";
my $uo = Bio::EnsEMBL::UnmappedObject->new (
	-unmapped_object_id 	=> 0,
       	-unmapped_reason_id 	=> 0,
	-type           	=> $test_type,
	-analysis       	=> $analysis,
	-external_db_id 	=> $test_ex_db_id,
	-identifier     	=> $test_identifier,
	-query_score		=> $test_query_score,
	-target_score		=> $test_target_score,
	-ensembl_id		=> $test_ensembl_id,
	-ensembl_object_type 	=> $test_ensembl_object_type,
	-summary		=> $test_summary,
	-full_desc      	=> $test_desc);
	 

$uma->store($uo);

my @objects = @{$uma->fetch_all()};
ok(scalar(@objects) == 5);

$uo = @{$uma->fetch_by_identifier($test_identifier)}[0];

ok($test_identifier eq $uo->identifier());
ok($test_ex_db_id == $uo->external_db_id());
ok($test_query_score == $uo->query_score());
ok($test_target_score == $uo->target_score());
ok($test_ensembl_id == $uo->ensembl_id());
ok($test_ensembl_object_type eq $uo->ensembl_object_type());
ok($test_summary eq $uo->summary());
ok($test_desc eq $uo->description());
ok($test_type eq $uo->type());

$uo->identifier('X1200');
ok($test_identifier ne $uo->identifier());
$uo->external_db_id(4200);
ok($test_ex_db_id != $uo->external_db_id());
$uo->query_score(45.3);
ok($test_query_score != $uo->query_score());
$uo->target_score(10.7);
ok($test_target_score != $uo->target_score());
$uo->summary("short fake one");
ok($test_summary ne $uo->summary());
$uo->description('long fake description');
ok($test_desc ne $uo->description());
$uo->ensembl_id(98765);
ok($test_ensembl_id != $uo->ensembl_id());
$uo->ensembl_object_type("Transcript");
ok($test_ensembl_object_type ne $uo->ensembl_object_type());


ok($uo->unmapped_reason_id());
ok($uo->dbID());

# now try minumum input test
$test_identifier ="X4444";
$uo = Bio::EnsEMBL::UnmappedObject->new (
	-type           => $test_type,
	-analysis       => $analysis,
	-external_db_id => $test_ex_db_id,
	-identifier     => $test_identifier,
	-summary	=> $test_summary,
	-full_desc      => $test_desc);


$uma->store($uo);

$uo = @{$uma->fetch_by_identifier($test_identifier)}[0];

ok($test_identifier eq $uo->identifier());
ok($test_ex_db_id == $uo->external_db_id());
ok($test_summary eq $uo->summary());
ok($test_desc eq $uo->description());
ok($test_type eq $uo->type());

ok($uo->unmapped_reason_id());
ok($uo->dbID());



$multi->restore('core', 'unmapped_object');
$multi->restore('core', 'unmapped_reason');


done_testing();
