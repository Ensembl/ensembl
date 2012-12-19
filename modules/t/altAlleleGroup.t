use strict;
use warnings;

use Bio::EnsEMBL::AltAlleleGroup;
use Bio::EnsEMBL::DBSQL::AltAlleleGroupAdaptor;

our $verbose = 0; #set to 1 to turn on debug printouts
use Test::More;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;

# Tests for basic methods on fake data

my $group = Bio::EnsEMBL::AltAlleleGroup->new(
    -MEMBERS => [
        [qw( id1 0 MANUAL )],
        [qw( id2 1 MANUAL )],
        [qw( id3 0 MANUAL )],
    ]
);

my $other_group = Bio::EnsEMBL::AltAlleleGroup->new(
    -MEMBERS => [
        [qw( id4 0 AUTOMATIC )],
        [qw( id5 0 AUTOMATIC )],
        [qw( id6 0 AUTOMATIC )],
    ]
);

my $id_list = $group->get_all_Gene_ids;

is_deeply($id_list,[qw( id1 id2 id3 )], "Get all Gene IDs within a prefab group");

my $id = $group->ref_Gene_id();
is($id,'id2',"Get reference Gene ID");
is($other_group->ref_Gene_id(),undef,"Test behaviour without a reference Gene set");

# Tests for the adaptor

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
ok(1);

my $db = $multi->get_DBAdaptor( 'core' );
ok($db);

my $aaga = $db->get_adaptor('human','core','AltAlleleGroup');

my $group_list = $aaga->fetch_all_Groups;
diag $group_list;

done_testing();