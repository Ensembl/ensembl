

## We start with some black magic to print on failure.
BEGIN { $| = 1; print "1..3\n";
	use vars qw($loaded); }

END {print "not ok 1\n" unless $loaded;}

use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::EMBLLOAD::Obj;

use Bio::SeqIO;


use lib 't';
use EnsTestDB;
$loaded = 1;
print "ok 1\n";    # 1st test passes.
    
my $ens_test = EnsTestDB->new();
    
# Load some data into the db
$ens_test->do_sql_file("t/staticgoldenpath.dump");

# Get an EnsEMBL db object for the test db
my $db = $ens_test->get_DBSQL_Obj;

#$Bio::EnsEMBL::FeatureFactory::USE_PERL_ONLY = 1;

print "ok 2\n";    
$db->static_golden_path_type('UCSC');

$stadaptor = $db->get_StaticGoldenPathAdaptor();

$vc2 = $stadaptor->fetch_VirtualContig_by_chr_name('chr2');

$seqout = Bio::SeqIO->new(-file => '>t/vc2embl.out' , '-format' => 'embl' );

$seqout->write_seq($vc2);

print "ok 3\n";

#unlink("t/vc2embl.out");
