use strict;
use Bio::EnsEMBL::Test::TestUtils;

use Bio::EnsEMBL::Test::MultiTestDB;


BEGIN { $| = 1;
	use Test;
	plan tests => 14;
}

my $multi_db = Bio::EnsEMBL::Test::MultiTestDB->new;
my $db = $multi_db->get_DBAdaptor('core');

my $verbose = 0;

# Test Creation

my $rfa = $db->get_RegulatorySearchRegionAdaptor();

ok(ref($rfa) && $rfa->isa('Bio::EnsEMBL::DBSQL::RegulatorySearchRegionAdaptor'));

#
## retrieve a search region via dbID
#
debug('---- fetch_by_dbID (default coords) ----');
my $feat = $rfa->fetch_by_dbID(2);
ok($feat->dbID == 2);
ok($feat->slice->seq_region_name() eq 'X');
ok($feat->start == 99639860);
ok($feat->end == 99646074);
ok($feat->strand == 1);

# List_dbids
my $ids = $rfa->list_dbIDs();
ok (@{$ids});

# test fetch_by_name
my $rsr = $rfa->fetch_by_name("CisRed_Search_11");
ok($rsr->dbID() == 2);

# Test store

$multi_db->save('core', 'regulatory_search_region');

my $slice = $db->get_SliceAdaptor()->fetch_by_seq_region_id(469293);

my $analysis = $db->get_AnalysisAdaptor->fetch_by_logic_name('cisred_search');

my $rf = Bio::EnsEMBL::RegulatorySearchRegion->new(-name      => 'new_search_region',
					      -start     => 100,
					      -end       => 220,
					      -strand    => -1,
					      -slice     => $slice,
					      -analysis  => $analysis,
					      -ensembl_object_type => "Gene",
					      -ensembl_object_id   => 18258,
					      -adaptor   => $rfa);

$rfa->store($rf);

ok($rf->dbID && $rf->adaptor());

$rf = $rfa->fetch_by_dbID($rf->dbID());

ok($rf->name()  eq 'new_search_region');
ok($rf->start() eq 100);
ok($rf->analysis()->logic_name() eq 'cisred_search');

$multi_db->restore('core', 'regulatory_search_region');
