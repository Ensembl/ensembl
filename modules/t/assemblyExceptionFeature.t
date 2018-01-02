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

use Test::More;
use Test::Warnings;

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::AssemblyExceptionFeature;

our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new;

# get a core DBAdaptor
#
my $dba = $multi->get_DBAdaptor("patch");
my $aefa = $dba->get_AssemblyExceptionFeatureAdaptor();
ok($aefa);

#
# 1 create a new AssemblyExceptionFeature
#
my $aef = new Bio::EnsEMBL::AssemblyExceptionFeature;
ok($aef);

#
# test the basic getter and setters
#

# start
ok(test_getter_setter($aef,'start',10));

# end
ok(test_getter_setter($aef,'end',14));

# type
ok(test_getter_setter($aef,'type', 'HAP'));

# check adaptor attaching
$aef->adaptor($aefa);
is(ref($aef->adaptor), 'Bio::EnsEMBL::DBSQL::AssemblyExceptionFeatureAdaptor', "Created AssemblyExceptionFeature adaptor");

# fetch all
my $chr_slice = $dba->get_SliceAdaptor->fetch_by_region('chromosome', 
                                                        'Y');
my $ref_slice = $dba->get_SliceAdaptor->fetch_by_region('chromosome', '6');
my $patch_slice = $dba->get_SliceAdaptor->fetch_by_region('supercontig', 'HG1304_PATCH');

my @features = @{$aefa->fetch_all_by_Slice($chr_slice)};

is(@features, 2, "Fetched one assembly exception feature for Y");

my @ref_features = @{$aefa->fetch_all_by_Slice($ref_slice)};
is(@ref_features, 1, "Fetched one assembly exception features for chromosome 6");

my @patch_features = @{ $aefa->fetch_all_by_Slice($patch_slice) };
is(@patch_features, 1, "Fetched one assembly exception for HG1304_PATCH");

foreach my $f (@features) {
  debug( "Feature: " . $f->slice->seq_region_name . " " . 
         $f->start . " " . $f->end . " " . $f->type);
  my $as = $f->alternate_slice();
  debug(" Alternate slice: " . $as->seq_region_name . " " . 
        $as->start . " " . $as->end);
}

foreach my $f (@ref_features) {
  debug( "Feature: " . $f->slice->seq_region_name . " " .
         $f->start . " " . $f->end . " " . $f->type);
  my $as = $f->alternate_slice();
  debug(" Alternate slice: " . $as->seq_region_name . " " .
        $as->start . " " . $as->end);
}

foreach my $f (@patch_features) {
  debug( "Feature: " . $f->slice->seq_region_name . " " .
         $f->start . " " . $f->end . " " . $f->type);
  my $as = $f->alternate_slice();
  debug(" Alternate slice: " . $as->seq_region_name . " " .
        $as->start . " " . $as->end);
}

my ($f) = @features;
is($f->display_id, $f->alternate_slice->seq_region_name, "Feature display id matches feature's alternate slice name");


my $feat = $aefa->fetch_by_dbID(1);

is($feat->dbID(), 1, "Feature dbID is 1");

#check we can store assembly exception features
my $aef_store = new Bio::EnsEMBL::AssemblyExceptionFeature();
my $aef_store2 = new Bio::EnsEMBL::AssemblyExceptionFeature();
#get ref slice
my $ref_slice = $dba->get_SliceAdaptor->fetch_by_region('chromosome', 6);
#prepare first object, the haplotype
$aef_store->start(1500);
$aef_store->end(35000);
$aef_store->type('HAP');
$aef_store->slice($chr_slice);
$aef_store->alternate_slice($ref_slice);
#prepare second object, the ref region to be substituted
$aef_store2->start(4500);
$aef_store2->end(38000);
$aef_store2->type('HAP REF');
$aef_store2->slice($ref_slice);
$aef_store2->alternate_slice($chr_slice);

my $asx_id = $aefa->store($aef_store,$aef_store2);

my $aef_new = $aefa->fetch_by_dbID($asx_id);

is($aef_new->dbID, $asx_id, "Assembly exception new dbID matches the stored id");

$aefa->remove($aef_store);

done_testing();
