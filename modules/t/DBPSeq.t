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
BEGIN { $| = 1; print "1..9\n"; 
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
$ens_test->do_sql_file("t/DBPseq.dump");
    
# Get an EnsEMBL db object for the test db
my $db = $ens_test->get_DBSQL_Obj;
print "ok 2\n";    

#Get Contig
my $contig=$db->get_Contig('AC021078.00006');
print "ok 3\n";

#Get DBPrimarySeq for this contig
my $dbseq = $contig->primary_seq();
print "ok 4\n";

#Check that DBPrimarySeq methods work
my $contig_id = $dbseq->contig_internal_id;

print "ok 5\n";
#print STDERR "Internal id of $conf{'contig'} is $contig_id\n";

my $dna_id = $dbseq->dna_id;
my $primary_id = $dbseq->primary_id;
my $display_id = $dbseq->display_id;
my $id = $dbseq->id;
print "ok 6\n";
#print STDERR "Dna id of $conf{'contig'} is $dna_id\n";

my $seq = $dbseq->seq();
print "ok 7\n";


my $subseq = $dbseq->subseq(2,6);
if( $subseq ne 'CCGAT' ) {
  print "not ok 8\n";
} else {
  print "ok 8\n";
}

my $length = $dbseq->length();
if( $length != 1222 ) {
   print "not ok 9\n";	    
} else {
   print "ok 9\n";	    
}





				      
