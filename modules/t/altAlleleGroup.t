use strict;
use warnings;

use Bio::EnsEMBL::AltAlleleGroup;
use Bio::EnsEMBL::DBSQL::AltAlleleGroupAdaptor;

use Test::More;
use Test::Exception;
use Test::Differences;
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

# Tests for methods applied to test db

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
ok(1);
my $db = $multi->get_DBAdaptor( 'core' );
ok($db);

my $aaga = $db->get_adaptor('AltAlleleGroup');
# Test data consists of a single group, of type AUTOMATIC, with a reference Allele and 3 others

my $group_list = $aaga->fetch_all_Groups;
my $aag = $group_list->[0];
is($aag->rep_Gene_id, 18256,"Check for correct selection of allele");

is_deeply ($aag->get_all_Gene_ids,[18256,18257,18258,18259],"Check group members");
is_deeply ($aag->get_all_Gene_ids('no ref'),[18257,18258,18259],"Test effect of excluding reference gene");

my $gene = $aag->get_representative_Gene;
is($gene->stable_id, 'ENSG00000131044',"Ensure both correct instantiation of Gene and ID thereof");

my $gene_list = $aag->get_all_Genes;
$gene = $gene_list->[0];
is (ref($gene),'Bio::EnsEMBL::Gene',"Returned object type from get_all_Genes");
ok ($gene->dbID == 18256,"Ensure Gene objects acquire correct information");

# test fetch_all_Groups_by_type
$group_list = $aaga->fetch_all_Groups_by_type('UNLIKELY STRING');
ok(scalar(@{$group_list}) == 0,"Try outlandish typed group lookup");
$group_list = $aaga->fetch_all_Groups_by_type('HAS_CODING_POTENTIAL');
ok(scalar(@$group_list) == 1,"Try known group type lookup");

# fetch_Group_by_id
$aag = $group_list->[0];
my $group_id = $aag->dbID;
my $new_aag = $aaga->fetch_Group_by_id($group_id);

is_deeply($aag,$new_aag,"Compare previously fetched group with group found by using the dbID");

$aag = $aaga->fetch_Group_by_id(undef);

ok(!defined($aag),"See what happens if no argument is given");

# fetch_Group_by_Gene_dbID
$aag = $aaga->fetch_all_Groups->[0];
$new_aag = $aaga->fetch_Group_by_Gene_dbID(18257);
is_deeply($new_aag,$aag,"Check single gene ID returns entire group correctly");

# check store method
my $dbID = $aaga->store($group);
ok($dbID);

my $aag2 = $aaga->fetch_Group_by_id($dbID);
is_deeply($aag2->get_all_members,$group->get_all_members,"Compare stored with original");
ok( $aaga->remove($aag2) );
$group->dbID($dbID);
$group->rep_Gene_id(1);
note("dbID = ".$group->dbID());
my $new_dbID = $aaga->update($group);
note($new_dbID);
$aag2 = $aaga->fetch_Group_by_id($dbID);

is_deeply($aag2->get_all_Gene_ids,[1,2,3], "Update and re-retrieve the same AltAlleleGroup") or diag explain $aag2->get_all_Gene_ids;
eq_or_diff($aag2->get_all_Gene_ids,[1,2,3], "Update and re-retrieve the same AltAlleleGroup");

# Vaguely verify the AltAlleleGroupAdaptor's fetch_all_Groups with a multispecies database
# Proper test data is hard to fabricate and no samples exist.

$aaga->db->is_multispecies(1);
$group_list = $aaga->fetch_all_Groups;
$aag = $group_list->[0];

ok(scalar(@$group_list) == 1, "Pretend multi-species fetch returns same groups as normal.");


done_testing();
