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
BEGIN { $| = 1; print "1..4\n"; 
	use vars qw($loaded); }
END {print "not ok 1\n" unless $loaded;}


use Bio::EnsEMBL::DBSQL::SymmetricContigFeatureContainer;

use lib 't';
use EnsTestDB;
$loaded = 1;
print "ok 1\n";    # 1st test passes.
    
my $ens_test = EnsTestDB->new();
$ens_test->do_sql_file('t/symmetriccf.dump');
my $db = $ens_test->get_DBSQL_Obj;

my $symadp = Bio::EnsEMBL::DBSQL::SymmetricContigFeatureContainer->new($db);

print "ok 2\n";



$fp = Bio::EnsEMBL::FeatureFactory->new_feature_pair();

$fp->set_all_fields(10,20,1,400,1,'symmetric',1,
		     40,50,1,400,2,'symmetric',2);

$symadp->write_FeaturePair_List($fp);

print "ok 3\n";

($fp) = $symadp->get_FeaturePair_list_by_rawcontig_id(1);

if( $fp->start == 10 && $fp->end == 20 && $fp->strand == 1 && $fp->hstart == 40 && $fp->hend == 50 && $fp->score == 400) {
      print "ok 4\n";
} else {
       print " not ok 4\n";
}



