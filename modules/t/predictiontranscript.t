use lib 't';

BEGIN { $| = 1;  
	use Test;
	plan tests => 10;
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
my $analysis = $db->get_AnalysisAdaptor->fetch_by_newest_logic_name("dummy-genscan");

ok($analysis);
ok($contig);

$prediction_f_ad = $db->get_PredictionTranscriptAdaptor();

my $ptrans = Bio::EnsEMBL::PredictionTranscript->new();

my $exon1 = Bio::EnsEMBL::Exon->new();
my $exon2 = Bio::EnsEMBL::Exon->new();

$exon1->start(19);
$exon1->end(189);
$exon1->strand(-1);
$exon1->phase(1);
$exon1->attach_seq($contig->primary_seq);
$exon1->contig($contig);
$exon2->start(269);
$exon2->end(457);
$exon2->strand(-1);
$exon2->phase(2);
$exon2->attach_seq($contig->primary_seq);
$exon2->contig($contig);
$exon1->analysis($analysis);
$exon2->analysis($analysis);
$ptrans->add_Exon($exon1);
$ptrans->add_Exon($exon2);
$ptrans->analysis($analysis);
@exons = $ptrans->get_all_Exons;



ok($ptrans);
@exons = $ptrans->get_all_Exons;


$prediction_f_ad->store($ptrans);
@exons = $ptrans->get_all_Exons;


ok(1);

my $out_preds = $prediction_f_ad->fetch_by_Contig($contig);

my $out_pred = $out_preds->[0];
ok($out_pred);

@exons = $out_pred->get_all_Exons;

ok(@exons);
ok($exons[0]->phase == 1);
ok($exons[1]->phase == 2);
