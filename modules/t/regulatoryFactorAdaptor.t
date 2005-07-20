use strict;
use Bio::EnsEMBL::Test::TestUtils;

use Bio::EnsEMBL::Test::MultiTestDB;


BEGIN { $| = 1;
	use Test;
	plan tests => 11;
}

my $multi_db = Bio::EnsEMBL::Test::MultiTestDB->new;
my $db = $multi_db->get_DBAdaptor('core');

my $verbose = 0;

# Test Creation

my $rma = $db->get_RegulatoryFactorAdaptor();

ok(ref($rma) && $rma->isa('Bio::EnsEMBL::DBSQL::RegulatoryFactorAdaptor'));

#
# Test fetch_by_dbID
#

my $rm = $rma->fetch_by_dbID(1);
ok($rm->name() eq 'hsa-miR-23b');
ok($rm->dbID == 1);
ok($rm->type eq 'miRNA_target');

#
# Test fetch_by_name
#
$rm = $rma->fetch_by_name('hsa-miR-23b');
ok($rm->name() eq 'hsa-miR-23b');
ok($rm->dbID() == 1);
ok($rm->type() eq 'miRNA_target');

#
# Test fetch_all_by_type
#
ok(@{$rma->fetch_all_by_type('transcription_factor')} == 2);

#
# Test store
#

$multi_db->save('core', 'regulatory_factor');

$rm = Bio::EnsEMBL::RegulatoryFactor->new(-NAME => 'test_store',
					 -TYPE => 'transcription_factor');

$rma->store($rm);

ok($rm->dbID && $rm->adaptor());

$rm = $rma->fetch_by_dbID($rm->dbID);

ok($rm->name eq 'test_store');
ok($rm->type eq 'transcription_factor');

$multi_db->restore('core', 'regulatory_factor');
