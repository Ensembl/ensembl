use lib 't';

BEGIN { $| = 1;  
	use Test;
	plan tests => 11;
}

my $loaded = 0;
END {print "not ok 1\n" unless $loaded;}

use EnsTestDB;
use Bio::EnsEMBL::DBLoader;



$loaded = 1;

ok(1);

# Database will be dropped when this
# object goes out of scope
my $ens_test = EnsTestDB->new;

$ens_test->do_sql_file("t/minidatabase.dump");

ok($ens_test);



my $db = $ens_test->get_DBSQL_Obj;
$cadp = $db->get_RawContigAdaptor();
$contig = $cadp->fetch_by_dbID(1);
my $analysis = $db->get_AnalysisAdaptor->fetch_by_logic_name("dummy-repeatmask");
print STDERR "ANALYSIS ".$analysis."\n";
ok($analysis);
ok($contig);

$repeat_f_ad = $db->get_RepeatFeatureAdaptor();
$repeat_c_ad = $db->get_RepeatConsensusAdaptor();


print STDERR "DBID ".$analysis->dbID."\n";

my $repeat_consensus = Bio::EnsEMBL::RepeatConsensus->new();

$repeat_consensus->length(10);
$repeat_consensus->repeat_class('dummy');
$repeat_consensus->name('dummy');
$repeat_consensus->repeat_consensus('ATGCATGCAT');

ok($repeat_consensus);

$repeat_c_ad->store($repeat_consensus);

ok(1);

my $repeat_feature = Bio::EnsEMBL::RepeatFeature->new();

$repeat_feature->contig_id($contig->dbID);
$repeat_feature->start(26);
$repeat_feature->end(65);
$repeat_feature->strand(1);
$repeat_feature->hstart(6);
$repeat_feature->hend(45);
$repeat_feature->score(100);
$repeat_feature->analysis($analysis);
$repeat_feature->repeat_consensus($repeat_consensus);

ok($repeat_feature);

$repeat_f_ad->store($contig->dbID, $repeat_feature);

ok(1);


my @repeats = $repeat_f_ad->fetch_by_RawContig($contig);

my $repeat = $repeats[0];

ok($repeat);

ok($repeat->start == 26);
ok($repeat->hend == 45);

