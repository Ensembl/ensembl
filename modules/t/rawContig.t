
use lib 't';

BEGIN { $| = 1;  
	use Test;
	plan tests => 8;
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

ok($contig);


$contig = $cadp->fetch_by_name('dummy-contig-1');

ok($contig->dbID == 1);

ok($contig->name eq 'dummy-contig-1');

$clone = $contig->clone();

ok($clone);

ok($clone->id eq 'dummy-embl-acc');

ok($contig->seq);
