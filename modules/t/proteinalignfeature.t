
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

ok($ens_test);

my @feats;

my $db = $ens_test->get_DBSQL_Obj;

$pep_f_ad = $db->get_ProteinAlignFeatureAdaptor();

my $feature1 = new Bio::EnsEMBL::SeqFeature();
$feature1->start(50);
$feature1->end(115);
$feature1->strand(-1);
$feature1->score(41);
$feature1->seqname(1);
$feature1->percent_id(100);

my $feature2 = new Bio::EnsEMBL::SeqFeature();
$feature2->start(5);
$feature2->end(26);
$feature2->seqname('dummy-hid');
$feature2->score(41);
$feature2->strand(1);
$feature2->percent_id(100);

$fp = new Bio::EnsEMBL::FeaturePair(-feature1 => $feature1,
				    -feature2 => $feature2);



push(@feats, $fp);


$feature1 = new Bio::EnsEMBL::SeqFeature();
$feature1->start(07);
$feature1->end(39);
$feature1->strand(-1);
$feature1->score(41);
$feature1->seqname(1);
$feature1->percent_id(100);

$feature2 = new Bio::EnsEMBL::SeqFeature();
$feature2->start(27);
$feature2->end(37);
$feature2->seqname('dummy-hid');
$feature2->score(41);
$feature2->strand(1);
$feature2->percent_id(100);

$fp2 = new Bio::EnsEMBL::FeaturePair(-feature1 => $feature1,		
				     -feature2 => $feature2);			

	
push(@feats, $fp2);

$pepf = Bio::EnsEMBL::DnaPepAlignFeature->new( -features => \@feats );

$pepf->analysis($db->get_AnalysisAdaptor->fetch_by_logic_name("dummy-blast"));

ok($pepf);

$pep_f_ad->store(1, $pepf);

ok(1);

($outf) = $pep_f_ad->fetch_by_contig_id(1);

ok($outf);  


ok($outf->start == 07);
ok($outf->end == 115);

