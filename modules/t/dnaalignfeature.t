

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

$dna_f_ad = $db->get_DnaAlignFeatureAdaptor();




$feature1 = new Bio::EnsEMBL::SeqFeature();
$feature1->start(5);
$feature1->end  (7);
$feature1->strand(1);
$feature1->score(10);
$feature1->seqname(1);
#$feature1->analysis($self->analysis);

$feature2 = new Bio::EnsEMBL::SeqFeature();
$feature2->start  (105);
$feature2->end    (107);
$feature2->strand (1);
$feature2->score  (10);
$feature2->seqname('dummy-hid');
$fp = new Bio::EnsEMBL::FeaturePair(-feature1 => $feature1,
				    -feature2 => $feature2);

push(@feats,$fp);


$feature1 = new Bio::EnsEMBL::SeqFeature();
$feature1->start(10);
$feature1->end  (14);
$feature1->strand(1);
$feature1->score(10);
$feature1->seqname(1);
#$feature1->analysis($self->analysis);

$feature2 = new Bio::EnsEMBL::SeqFeature();
$feature2->start  (106);
$feature2->end    (110);
$feature2->strand (1);
$feature2->score  (10);
$feature2->seqname('dummy-hid');
$fp2 = new Bio::EnsEMBL::FeaturePair(-feature1 => $feature1,
				    -feature2 => $feature2);


push(@feats,$fp2);


$dnaf = Bio::EnsEMBL::DnaDnaAlignFeature->new( -features => \@feats );


$dnaf->seqname(1);
$dnaf->hseqname('dummy-hid');

$dnaf->analysis($db->get_AnalysisAdaptor->fetch_by_logic_name("dummy-blast"));

ok($dnaf);

$dna_f_ad->store(1,$dnaf);

ok(1);

($outf) = $dna_f_ad->fetch_by_contig_id(1);

ok($outf);
ok($outf->start == 5);
ok($outf->end   == 14);

ok( scalar($outf->ungapped_features) == 2);

($outf) = $dna_f_ad->fetch_by_assembly_location(1,400,1,'NCBI_28');

ok($outf);
ok( scalar($outf->ungapped_features) == 2);
