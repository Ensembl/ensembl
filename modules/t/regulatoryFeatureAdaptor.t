use strict;
use Bio::EnsEMBL::Test::TestUtils;

use Bio::EnsEMBL::Test::MultiTestDB;


BEGIN { $| = 1;
	use Test;
	plan tests => 13;
}

my $multi_db = Bio::EnsEMBL::Test::MultiTestDB->new;
my $db = $multi_db->get_DBAdaptor('core');

my $verbose = 0;

# Test Creation

my $rfa = $db->get_RegulatoryFeatureAdaptor();

ok(ref($rfa) && $rfa->isa('Bio::EnsEMBL::DBSQL::RegulatoryFeatureAdaptor'));

#
## retrieve a feature via dbID
#
debug('---- fetch_by_dbID (default coords) ----');
my $feat = $rfa->fetch_by_dbID(2);
ok($feat->dbID == 2);
ok($feat->slice->seq_region_name() eq '3');
ok($feat->start == 49546683);
ok($feat->end == 49546702);
ok($feat->strand == 1);

# test fetch all by motif
my $rma = $db->get_RegulatoryMotifAdaptor();
my $motif = $rma->fetch_by_name('motif1');
ok($motif);
ok(@{$rfa->fetch_all_by_motif($motif)} == 10);

# List_dbidx
my $ids = $rfa->list_dbIDs();
ok (@{$ids});

# test fetch_all_by_transcript
my $transcript =  $db->get_TranscriptAdaptor()->fetch_by_dbID(21740);
ok(@{$rfa->fetch_all_by_transcript($transcript)} == 13);

# Test store

$multi_db->save('core', 'regulatory_feature');

my $slice = $db->get_SliceAdaptor()->fetch_by_seq_region_id(469286);

my $analysis = $db->get_AnalysisAdaptor->fetch_by_logic_name('miRanda');

my $rm = Bio::EnsEMBL::RegulatoryMotif->new(-name => 'motif1',
					    -type => 'miRNA_target');

my $rf = Bio::EnsEMBL::RegulatoryFeature->new(-name      => 'hsa-miR-108',
					      -start     => 100,
					      -end       => 220,
					      -strand    => -1,
					      -slice     => $slice,
					      -analysis  => $analysis,
					      -motif     => $motif,
					      -influence => 'positive',
					      -adaptor   => $rfa);

$rfa->store($rf);

ok($rf->dbID && $rf->adaptor());

$rf = $rfa->fetch_by_dbID($rf->dbID());

ok($rf->name()  eq 'hsa-miR-108');
ok($rf->start() eq 100);
ok($rf->analysis()->logic_name() eq 'miRanda');

$multi_db->restore('core', 'regulatory_feature');
