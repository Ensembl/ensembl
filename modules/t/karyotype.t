## We start with some black magic to print on failure.
BEGIN { $| = 1; print "1..3\n";
	use vars qw($loaded); }

END {print "not ok 1\n" unless $loaded;}

use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::DBSQL::KaryotypeAdaptor;

use lib 't';
use EnsTestDB;
$loaded = 1;
print "ok 1\n";    # 1st test passes.

my $ens_test = EnsTestDB->new();
    
# Load some data into the db
$ens_test->do_sql_file("t/karyotype.dump");


# Get an EnsEMBL db object for the test db
my $db = $ens_test->get_DBSQL_Obj;
print "ok 2\n";    

$kadp = Bio::EnsEMBL::DBSQL::KaryotypeAdaptor->new($db);

$band = $kadp->get_band_label_by_position('chr1',1000);

if( $band eq 'p36.33') {
	print "ok 3\n";
} else {
	print "not ok 3\n";
}
