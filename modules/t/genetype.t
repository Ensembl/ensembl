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
BEGIN { $| = 1; print "1..5\n"; 
	use vars qw($loaded); }

END {print "not ok 1\n" unless $loaded;}

use lib 't';
use EnsTestDB;
$loaded = 1;
print "ok 1\n";    # 1st test passes.
    
my $ens_test = EnsTestDB->new();
    
# Load some data into the db
$ens_test->do_sql_file("t/genetype.dump");
    
# Get an EnsEMBL db object for the test db
my $db = $ens_test->get_DBSQL_Obj;
print "ok 2\n";    



my $gene_obj = $db->gene_Obj;

my $gene=$gene_obj->get('ENSG00000003941');


if ($gene){print "ok 3\n";}
else {print "Not ok 3\n";}



my $contig=$db->get_Contig('AB000381.00001');

my @genes=$contig->get_Genes_by_Type('ensembl');
if (scalar @genes !=0){print "ok 4\n";}
else {print "Not ok 4\n";}


my $vc = Bio::EnsEMBL::DB::VirtualContig->new( -focuscontig => $contig,
                                              -focusposition => 10000,
                                              -ori => 1,
                                              -left => 2000000,
                                              -right => 2000000
                                              );


my @genes=$vc->get_Genes_by_Type('ensembl');
if (scalar @genes !=0){print "ok 5\n";}
else {print "Not ok 5\n";}


