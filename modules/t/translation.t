use lib 't';
use TestUtils qw(test_getter_setter);
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Exon;

BEGIN { $| = 1;  
	use Test;
	plan tests => 10;
}

my $loaded = 0;
END {print "not ok 1\n" unless $loaded;}

my $verbose = 0;

use MultiTestDB;

my $multi = MultiTestDB->new();

$loaded = 1;

ok(1);

my $db = $multi->get_DBAdaptor( 'core' );

my $t = Bio::EnsEMBL::Translation->new();

ok($t);

ok(test_getter_setter($t,'stable_id',1));

ok(test_getter_setter($t,'dbID',3));

ok(test_getter_setter($t,'start',42));
ok(test_getter_setter($t,'end',50));

my $exon = Bio::EnsEMBL::Exon->new();
$exon->start(10);
$exon->end(20);
$exon->strand(1);
$exon->phase(0);
#$exon->contig( $contig );
$exon->end_phase( -1 );

$t->start_Exon($exon);
ok($t);

$t->end_Exon($exon);
ok($t);

my $ta = $db->get_TranslationAdaptor();
my $ids = $ta->list_dbIDs();
ok (@{$ids});

my $stable_ids = $ta->list_stable_ids();
ok (@{$stable_ids});

