use lib 't';

BEGIN { $| = 1;  
	use Test;
	plan tests => 7;
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
my $db = $ens_test->get_DBSQL_Obj;
ok($ens_test);

my $analysis_ad = $db->get_AnalysisAdaptor();

ok($analysis_ad);



my $analysis = Bio::EnsEMBL::Analysis->new();

$analysis->logic_name('dummy_analysis');
$analysis->db('dummy');
$analysis->program('dummy');
$analysis->gff_source('dummy');
$analysis->gff_feature('dummy');


ok($analysis);

$analysis_ad->store($analysis);


ok(1);


my $analysis_out = $analysis_ad->fetch_by_logic_name('dummy_analysis');


ok($analysis_out);
ok($analysis_out->db eq 'dummy');
