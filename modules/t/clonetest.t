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
BEGIN { $| = 1; print "1..19\n"; 
	use vars qw($loaded); }
END {print "not ok 1\n" unless $loaded;}


use lib 't';
use EnsTestDB;
use Bio::EnsEMBL::DBSQL::Clone;

$loaded = 1;
print "ok 1\n";    # 1st test passes.
    
my $ens_test = EnsTestDB->new();
    
# Load some data into the db
$ens_test->do_sql_file("t/clonetest.dump");
    
# Get an EnsEMBL db object for the test db
my $db = $ens_test->get_DBSQL_Obj;
print "ok 2\n";    

#Checking get clone  method exeption
my $clone;
eval {
    $clone=$db->get_Clone("wrong_id");
};
if ($@) {
    print "ok 3\n";
    
}
else {
    print "not ok 3\n";
    print STDERR "Trying to get a non existing clone not throw an exeption\n";
}

#Check if its get an existing clone
eval {
    $clone=$db->get_Clone("test3");
};
if ($@) {
    print "not ok 4\n";
    print STDERR "Could not get an existing clone using get_clone\n";
}
else {
    print "ok 4\n";
}

#Test if there is all of the features for clone
if ($clone->id eq "test3") {
    print "ok 5\n";   
}
else {
    print "not ok 5\n";
    print STDERR "clone->id does not give the expected value!\n";
}

if ($clone->dbID == 3) {
    print "ok 6\n";
    
}
else {
    print "not ok 6\n";
    print STDERR "clone->_internal_id does not give the expected value!\n";
}


if ($clone->embl_id eq "embl-test3" ) {
    print "ok 7\n";
}
else {
    print "not ok 7\n";
    print STDERR "clone->embl_id does not give the expected value!\n";

}

if ($clone->version == 1 ) {
    print "ok 8\n";
}
else {
    print "not ok 8\n";
    print STDERR "clone->version does not give the expected value!\n";
}

if ($clone->embl_version == 1 ) {
    print "ok 9\n";
}
else {
    print "not ok 9\n";
    print STDERR "clone->embl->version does not give the expected value!\n";
}

if ($clone->htg_phase == 1) {
    print "ok 10\n";
}
else {
    print "not ok 10\n";
    print STDERR "clone->htg_phase does not give the expected value!\n";
}

if ($clone->created == 971774887) {
    print "ok 11\n";
}
else {
    print "not ok 11\n";
    print STDERR "clone->created does not give the expected value!\n";
}

if ($clone->modified == 971774887) {
    print "ok 12\n";
}
else {
    print "not ok 12\n";
    print STDERR "clone->modified does not give the expected value!\n";
}




#Get all of the genes from the clone
my @genes =  $clone->get_all_Genes();

if ($genes[0]->id eq "gene_id3") {
    print "ok 13\n";
}
else {
    print "not ok 13\n";
    print STDERR "clone->get_all_Genes does not retrieve the expected gene!\n";
}

#Check the get contig method
@contig = $db->get_Contig(id_cont_test3);
if (scalar(@contig) != 1 || $contig[0]->id ne "id_cont_test3") {
    print "not ok 14\n";
    print STDERR "db->get_Contig does not retrieve the expected contig!\n";
}
else {
    print "ok 14\n";
    
}

 
@geneid = $clone->get_all_my_geneid();
if (scalar(@geneid) != 1 || $geneid[0] ne "gene_id3") {
    print "not ok 15\n";
    print STDERR "clone->get_all_my_geneid does not retrieve the expected gene id!\n";
}
else {
    print "ok 15\n";  
}

@allcontig = $clone->get_all_Contigs();
if (scalar(@allcontig) != 1 || $allcontig[0]->id ne "id_cont_test3") {
    print "not ok 16\n";
    print STDERR "clone->get_all_Contigs does not retrieve the expected contig!\n";
}
else {
    print "ok 16\n";  
}

@rowcontig = $clone->get_rawcontig_by_position(10);
if (scalar(@rowcontig) != 1 || $rowcontig[0]->id ne "id_cont_test3") {
    print "not ok 17\n";
    print STDERR "clone->get_rawcontig_by_position does not retrieve the expected contig\n";
}
else {
    print "ok 17\n";    
}

#Delete the current clone and test if this clone and the contigs have been deleted from the database
$clone->delete_by_dbID();

#Check if the clone has been deleted
eval {
    $clone=$db->get_Clone("test3");
};
if ($@) {
    print "ok 18\n";
}
else {
    print "not ok 18\n";
    print STDERR "clone->delete does not have deleted the clone!\n";
}

#Check if the contig has been deleted
eval {
@allcontig = $db->get_Contig(id_cont_test3);
};
if ($@) {
    print "ok 19\n";
}
else {
    print "not ok 19\n";
    print STDERR "clone delete does not have deleted the clone!\n";
}

