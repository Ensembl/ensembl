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
BEGIN { $| = 1; print "1..15\n"; 
	use vars qw($loaded); }

END {print "not ok 1\n" unless $loaded;}

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

my $gad =$db->get_GeneAdaptor;

#Checking get_all_Gene_id;

my @genes= $gad->list_stable_geneIds;


if ((scalar @genes == 1) && ($genes[0] eq 'ENSG1')) {
    print "ok 3\n";
}
else {
    print "not ok 3\n";
    print STDERR "Could not get gene ids using get_all_Gene_id\n";
}

#Checking get method exception
eval{
my $gene = $gad->fetch_by_dbID(10000);
};
if ($@ =~ /EXCEP/) {
    print "ok 4\n";
}
else {
    print "not ok 4\n";
    print STDERR "Trying to get a non-existing gene with gene_obj->get did not throw an exception! $@\n";
} 

#Checking get method (correct use)
my $gene = $gad->fetch_by_stable_id('ENSG1');
if ($gene->isa('Bio::EnsEMBL::Gene')) {
    print "ok 5\n";
}
else {
    print "not ok 5\n";
    print STDERR "Could not get the test gene from the database!\n";
}

#Check gene methods
if ($gene->stable_id eq 'ENSG1') {
    print "ok 6\n";
}
else {
    print "not ok 6\n";
    print "gene->stable_id does not give the expected value!\n";
}

if ($gene->version == 1) {
    print "ok 7\n";
}
else {
    print "not ok 7\n";
    print STDERR "gene->version does not give the expected value! - instead it give ",$gene->version,"\n";
}

# skipping date tests


    print "ok 8\n";
    print "ok 9\n";
    print "ok 10\n";

$expected = 'some description';
$desc = $gene->description;
if ( $desc eq $expected ) {
    print "ok 11\n";
}
else {
    print "not ok 11\n";
    print STDERR "expected $expected, got $desc\n";
} 


if ($gene->is_known == 1) {
    print "ok 12\n";
}
else {
    print "not ok 12\n";
    print STDERR "Gene contains DBLinks, but $gene->is_known does not return 1\n";
}
$ok=0;

$ok=0;
foreach my $exon ($gene->get_all_Exons) {
    if ($exon->stable_id =~ /exon-1|exon-2/) {
	print STDERR $exon->contig_id(),"\n",$exon->clone_id(),"\n";
	$ok++;
    }
}
if ($ok == 2) {
    print "ok 13\n";
}
else {
    print "not ok 13\n";
    print STDERR "$gene->get_all_Exons does not give expected exons $ok\n";
}


#Not testing add_transcript yet
$ok=0;
foreach my $trans ($gene->each_Transcript) {
    if ($trans->stable_id =~ /ENST1/) {
	$ok++;
    }
}
if ($ok == 1) {
    print "ok 14\n";
}
else {
    print "not ok 14\n";
}


# test removeing gene

$gad->remove( $gene );

if( ! defined $gene->dbID() ) {
	print "ok 15\n";
} else {
	print "not ok 15\n";
}

