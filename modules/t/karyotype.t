## We start with some black magic to print on failure.
BEGIN { $| = 1; print "1..3\n";
	use vars qw($loaded); }

END {print "not ok 1\n" unless $loaded;}

use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::DBSQL::KaryotypeBandAdaptor;

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

my $kadp = Bio::EnsEMBL::DBSQL::KaryotypeBandAdaptor->new($db);

my $bandobj = $kadp->fetch_by_chromosome_position('chr1',1000);

if ($bandobj){
	print "ok 3\n";
} else {
	print "not ok 3\n";
}

my $bandname = $bandobj->name;
if( $bandname eq 'p36.33') {
	print "ok 4\n";
} else {
	print "not ok 4\n";
}

my $bandchr = $bandobj->chromosome;
if( $bandchr eq 'chr1') {
	print "ok 5\n";
} else {
	print "not ok 5\n";
}

my $bandstart = $bandobj->start;
if( $bandstart == 0) {
	print "ok 6\n";
} else {
	print "not ok 6\n";
}

my $bandend = $bandobj->end;
if( $bandend == 2300000) {
	print "ok 7\n";
} else {
	print "not ok 7\n";
}

my $bandstain = $bandobj->stain;
if( $bandstain eq 'gneg') {
	print "ok 8\n";
} else {
	print "not ok 8\n";
}


$bandobj = $kadp->fetch_by_chromosome_name('chr1','p36.33');

if ($bandobj){
	print "ok 9\n";
} else {
	print "not ok 9\n";
}

my $bandname = $bandobj->name;
if( $bandname eq 'p36.33') {
	print "ok 10\n";
} else {
	print "not ok 10\n";
}

my $bandchr = $bandobj->chromosome;
if( $bandchr eq 'chr1') {
	print "ok 11\n";
} else {
	print "not ok 11\n";
}

my $bandstart = $bandobj->start;
if( $bandstart == 0) {
	print "ok 12\n";
} else {
	print "not ok 12\n";
}

my $bandend = $bandobj->end;
if( $bandend == 2300000) {
	print "ok 13\n";
} else {
	print "not ok 13\n";
}

my $bandstain = $bandobj->stain;
if( $bandstain eq 'gneg') {
	print "ok 14\n";
} else {
	print "not ok 14\n";
}


my @band = $kadp->fetch_all_by_chromosome('chr1');

if (scalar(@band)==2){
	print "ok 15\n";
} else {
	print "not ok 15\n";
}


