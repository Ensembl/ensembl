use lib 't';

BEGIN { $| = 1;  
	use Test;
	plan tests => 12;
}


use MultiTestDB;
use TestUtils qw(test_getter_setter debug);

our $verbose = 0;


my $multi = MultiTestDB->new();

ok(1);

my $db = $multi->get_DBAdaptor( 'core' );

$cadp = $db->get_RawContigAdaptor();

$contig = $cadp->fetch_by_dbID(319456);

my $analysis = $db->get_AnalysisAdaptor->fetch_by_logic_name("RepeatMask");

debug( "ANALYSIS ".$analysis );

ok($analysis);
ok($contig);

$repeat_f_ad = $db->get_RepeatFeatureAdaptor();
$repeat_c_ad = $db->get_RepeatConsensusAdaptor();


debug( "Analysis dbID ".$analysis->dbID );

my $repeat_consensus = Bio::EnsEMBL::RepeatConsensus->new();

$repeat_consensus->length(10);
$repeat_consensus->repeat_class('dummy');
$repeat_consensus->name('dummy');
$repeat_consensus->repeat_consensus('ATGCATGCAT');

ok($repeat_consensus);

$repeat_c_ad->store($repeat_consensus);

ok(1);

my $repeat_feature = Bio::EnsEMBL::RepeatFeature->new();

$repeat_feature->start(26);
$repeat_feature->end(65);
$repeat_feature->strand(1);
$repeat_feature->hstart(6);
$repeat_feature->hend(45);
$repeat_feature->score(100);
$repeat_feature->analysis($analysis);
$repeat_feature->repeat_consensus($repeat_consensus);
$repeat_feature->contig( $contig );

ok($repeat_feature);
$multi->hide( "core", "repeat_feature" );

$repeat_f_ad->store( $repeat_feature );


ok(1);

my $repeats = $repeat_f_ad->fetch_all_by_RawContig($contig);

my $repeat = $repeats->[0];

ok($repeat);

ok($repeat->start == 26);
ok($repeat->hend == 45);

my $dbID = $repeat->dbID;

my $r = $repeat_f_ad->fetch_by_dbID($dbID);

ok($r->dbID == $dbID && $r->start == 26 && $r->hend == 45);


# list_dbIDs
my $ids = $repeat_f_ad->list_dbIDs();
ok (@{$ids});

$multi->restore('core', 'repeat_feature');


