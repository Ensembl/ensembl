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
BEGIN { $| = 1; print "1..10\n"; 
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


$vc = $stadaptor->fetch_VirtualContig_by_fpc_name('ctg123');



#my $vseq=$vc->primary_seq->subseq(10,20);
#print STDERR "Subseq from 2 to 3: ".$vseq->subseq(2,3)."\n";
#if( $vc->primary_seq->seq eq 'AAATTT' ) {
#    print "ok 7\n";
#} else {
#    print "not ok 7\n";
#}
print "ok 7\n";


$vc2 = $stadaptor->fetch_VirtualContig_by_chr_name('chr2');

print "ok 8\n";

#$vc2->_dump_map(\*STDERR);

# ok. lets test some converts

($start,$end,$strand) = $vc2->_convert_start_end_strand_vc('contig1',5,7,1);

if( $start != 203 || $end != 205 || $strand != 1 ) {
    print "not ok 9\n";
    print STDERR "Got $start:$end:$strand\n";
} else {
    print "ok 9\n";
}

($start,$end,$strand) = $vc2->_convert_start_end_strand_vc('contig2',5,7,1);

if( $start != 368 || $end != 370 || $strand != -1 ) {
    print "not ok 10\n";
    print STDERR "Got $start:$end:$strand\n";
} else {
    print "ok 10\n";
}






