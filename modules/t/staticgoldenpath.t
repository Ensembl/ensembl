## Bioperl Test Harness Script for Modules
##
# CVS Version
# $Id$


# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

#-----------------------------------------------------------------------
## perl test harness expects the following output syntax only!
## 1..3
## ok 1  [not ok 1 (if test fails)]
## 2..3
## ok 2  [not ok 2 (if test fails)]
## 3..3
## ok 3  [not ok 3 (if test fails)]
##
## etc. etc. etc. (continue on for each tested function in the .t file)
#-----------------------------------------------------------------------


## We start with some black magic to print on failure.
BEGIN { $| = 1; print "1..8\n"; 
	use vars qw($loaded); }

END {print "not ok 1\n" unless $loaded;}

use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DBLoader;
use lib 't';
use EnsTestDB;
$loaded = 1;
print "ok 1\n";    # 1st test passes.
    
my $ens_test = EnsTestDB->new();
    
# Load some data into the db
$ens_test->do_sql_file("t/staticgoldenpath.dump");

# Get an EnsEMBL db object for the test db
my $db = $ens_test->get_DBSQL_Obj;
print "ok 2\n";    
$db->static_golden_path_type('UCSC');

$stadaptor = $db->get_StaticGoldenPathAdaptor();

@array = $stadaptor->fetch_RawContigs_by_fpc_name('ctg123');
if( scalar(@array) != 3 ) {
   print "not ok 3\n";
} else {
   print "ok 3\n";
}

$rc1 = shift @array;
if( $rc1->id ne 'contig1' ) {
   print "not ok 4\n";
} else {
   print "ok 4\n";
}

@array = $stadaptor->fetch_RawContigs_by_chr_name('chr2');
if( scalar(@array) != 3 ) {
   print "not ok 5\n";
} else {
   print "ok 5\n";
}

$rc1 = shift @array;
if( $rc1->id ne 'contig1' ) {
   print "not ok 6\n";
} else {
   print "ok 6\n";
}


$vc = $stadaptor->VirtualContig_by_fpc_name('ctg123');

print "ok 7\n";


$vc2 = $stadaptor->VirtualContig_by_chr_name('chr2');

print "ok 8\n";





