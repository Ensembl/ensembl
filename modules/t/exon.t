

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

# Exon specific tests

$exonad = $db->get_ExonAdaptor();

ok($exonad);

$exon = Bio::EnsEMBL::Exon->new();

$exon->start(10);
ok($exon->start == 10);

$exon->end(20);
ok($exon->end == 20);

$exon->strand(1);
ok($exon->strand == 1);

$exon->phase(0);
ok($exon->phase == 0);

$exon->contig_id(1);
ok($exon->contig_id == 1);

# should try to store (!)

$exonad->store($exon);

ok($exon->dbID);

# now test fetch_by_dbID

$newexon = $exonad->fetch_by_dbID($exon->dbID);

ok($newexon);




