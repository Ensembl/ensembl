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

use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DBLoader;

use lib 't';
use EnsTestDB;
$loaded = 1;
print "ok 1\n";    # 1st test passes.
    
my $ens_test = EnsTestDB->new();
    
# Load some data into the db
$ens_test->do_sql_file("t/geneget.dump");
    
# Get an EnsEMBL db object for the test db
my $db = $ens_test->get_DBSQL_Obj;
print "ok 2\n";    


$gene_obj = $db->gene_Obj();

eval {
    $gene = $gene_obj->get('ENSG1');
};

if ($@) {
    print "not ok 3\n";
}

else {
    print "ok 3\n";
}

@dblink = $gene->each_DBLink();

$dbl = shift @dblink;

if( !defined $dbl || $dbl->database ne 'swissprot' ) {
    print "not ok 4\n";
} else {
  print "ok 4\n";
}

print STDERR $gene->is_known,"\n";

#@trans = $gene->each_Transcript();
#$trans = shift @trans;






