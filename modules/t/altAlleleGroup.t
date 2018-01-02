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
use warnings;

use Bio::EnsEMBL::AltAlleleGroup;
use Bio::EnsEMBL::DBSQL::AltAlleleGroupAdaptor;

use Test::More;
use Test::Exception;
use Test::Differences;
use Test::Warnings;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;

# Tests for basic methods on fake data

my $group = Bio::EnsEMBL::AltAlleleGroup->new(
    -MEMBERS => [
        [ 1, {} ],
        [ 2, {IS_REPRESENTATIVE => 1}],
        [ 3, {}],
    ] # gene_id , type hash
);

my $other_group = Bio::EnsEMBL::AltAlleleGroup->new(
    -MEMBERS => [
        [4, {}],
        [5, {}],
        [6, {}],
    ]
);

ok($group->size == 3,"Size reporting correctly");

my $id_list = $group->get_all_Gene_ids;

is_deeply($id_list,[qw( 1 2 3 )], "Get all Gene IDs within a prefab group");

my $id = $group->rep_Gene_id();
is($id,'2',"Get reference Gene ID");
warns_like {
    is($other_group->rep_Gene_id(),undef,"Test behaviour without a reference Gene set");
} qr/No representative allele currently set/, 'Ensuring appropriate warning about lack of reference gene is emitted';

$group->unset_rep_Gene_id();
warns_like {
    is($group->rep_Gene_id(), undef,"Check successful unsetting.");
} qr/No representative allele currently set/, 'Ensuring appropriate warning about lack of reference gene is emitted';

throws_ok{$group->rep_Gene_id(5)} qr/Requested representative gene ID/, "Exception thrown when invalid ID given";
ok($group->rep_Gene_id(3) == 3,"Successful setting of reference gene");

$other_group->remove_all_members;
ok($other_group->size == 0, "Test remove_all_members");

# Test the methods which attempt to modify an AltAlleleGroup object
{
    my $copy_group = bless({%{$group}}, ref($group));
    ok($copy_group->contains_member(1), 'Group contains the member 1');
    ok(!$copy_group->contains_member(4), 'Group does not contain the member 4');
    is_deeply($copy_group->attribs(1), {}, 'Group member 1 attributes are empty');
    $copy_group->set_attribs(1, ['IN_CORRECTED_ASSEMBLY']);
    $copy_group->set_attribs(1, {'manually_assigned' => 1}); #lower case intentional
    $copy_group->set_attribs(1, 'IS_VALID_ALTERNATE');
    is_deeply($copy_group->attribs(1), { map { $_, 1} qw/IN_CORRECTED_ASSEMBLY IS_VALID_ALTERNATE MANUALLY_ASSIGNED/ }, 'Group member 1 attributes are empty');
    $copy_group->add_member(1,{});
    is_deeply($copy_group->attribs(1), { map { $_, 1} qw/IN_CORRECTED_ASSEMBLY IS_VALID_ALTERNATE MANUALLY_ASSIGNED/ }, 'Attributes were not replaced because the member was already in');
    $copy_group->remove_attribs(1, 'IS_VALID_ALTERNATE');
    $copy_group->remove_attribs(1, ['IN_CORRECTED_ASSEMBLY']);
    $copy_group->remove_attribs(1, {MANUALLY_ASSIGNED => 1});
    is_deeply($copy_group->attribs(1), {}, 'Group member 1 attributes are empty');

    $copy_group->remove_member(4);
    is($copy_group->size(), 3, 'Removing an non-existent ID does nothing');
    $copy_group->remove_member(1);
    is($copy_group->size(), 2, 'Removing ID 1 reduces our group size');
    ok(!$copy_group->contains_member(1), 'Group does not contain the member 1 post removal');
    $copy_group->add_member(1,{});
    is_deeply($copy_group->attribs(1), {}, 'Group member 1 attributes are empty post removal and addition');
}

# Tests for methods applied to test db

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db = $multi->get_DBAdaptor( 'core' );
ok($db, 'DB was retrieved');

my $aaga = $db->get_AltAlleleGroupAdaptor;
my $ga = $db->get_GeneAdaptor();
# Test data consists of a single group, of type AUTOMATIC, with a reference Allele and 3 others

my $group_list = $aaga->fetch_all;
my $aag = $group_list->[0];
is($aag->rep_Gene_id, 18256,"Check for correct selection of allele");

is_deeply ($aag->get_all_Gene_ids,[18256,18257,18258,18259],"Check group members");
is_deeply ($aag->get_all_Gene_ids('no ref'),[18257,18258,18259],"Test effect of excluding reference gene");

my $gene = $aag->get_representative_Gene;
is($gene->stable_id, 'ENSG00000131044',"Ensure both correct instantiation of Gene and ID thereof");

#Checking we can filter the members list by a gene object and an ID
is_deeply($aag->get_all_Gene_ids(undef, [18259]), [18256, 18257,18258], 'Filtering out a Gene by ID');
is_deeply($aag->get_all_Gene_ids(undef, [$ga->fetch_by_dbID(18259)]), [18256, 18257,18258], 'Filtering out a Gene by object');

my $gene_list = $aag->get_all_Genes;
$gene = $gene_list->[0];
is (ref($gene),'Bio::EnsEMBL::Gene',"Returned object type from get_all_Genes");
ok ($gene->dbID == 18256,"Ensure Gene objects acquire correct information");

# test fetch_all
$group_list = $aaga->fetch_all('UNLIKELY STRING');
ok(scalar(@{$group_list}) == 0,"Try outlandish typed group lookup");
$group_list = $aaga->fetch_all('HAS_CODING_POTENTIAL');
ok(scalar(@$group_list) == 1,"Try known group type lookup");

# fetch_Group_by_id
$aag = $group_list->[0];
my $group_id = $aag->dbID;
my $new_aag = $aaga->fetch_by_dbID($group_id);

is_deeply($aag, $new_aag, "Compare previously fetched group with group found by using the dbID");

$aag = $aaga->fetch_by_dbID(undef);

ok(!defined($aag), "See what happens if no argument is given");

# fetch_Group_by_Gene_dbID
$aag = $aaga->fetch_all->[0];
$new_aag = $aaga->fetch_by_gene_id(18257);
is_deeply($new_aag,$aag,"Check single gene ID returns entire group correctly");

# Saving the multi DBs
$multi->save('core', qw/alt_allele alt_allele_group alt_allele_attrib/);

# check store method
my $dbID = $aaga->store($group);
ok($dbID, 'A dbID was returned from the store method');

{
    my $aag2 = $aaga->fetch_by_dbID($dbID);
    is_deeply($aag2->get_all_members,$group->get_all_members,"Compare stored with original");
    $group->add_member(4, {});
    my $update_dbID = $aaga->update($group);
    cmp_ok($update_dbID, '==', $dbID, 'We should have kept the db id of the group');
}

$group->remove_member(4);
$aaga->remove($group);
ok(! defined $aaga->fetch_by_dbID($dbID), 'Using a deleted ID means no group returned');
my $new_dbID = $aaga->store($group);

cmp_ok($new_dbID, '!=', $dbID, 'Should have been assgined a new ID');
my $aag2 = $aaga->fetch_by_dbID($new_dbID);
my $gene_ids = $aag2->get_all_Gene_ids();
eq_or_diff($gene_ids,[1,2,3], "Update and re-retrieve the same AltAlleleGroup") or diag explain $gene_ids;

# Vaguely verify the AltAlleleGroupAdaptor's fetch_all_Groups with a multispecies database
# Proper test data is hard to fabricate and no samples exist.

$aaga->db->is_multispecies(1);
$group_list = $aaga->fetch_all;
$aag = $group_list->[0];

ok(scalar(@$group_list) == 1, "Pretend multi-species fetch returns same groups as normal.");


done_testing();
