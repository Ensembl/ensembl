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
BEGIN { $| = 1; print "1..7\n"; 
	use vars qw($loaded); }

END {print "not ok 1\n" unless $loaded;}

use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Exon;

$loaded = 1;
print "ok 1\n";    # 1st test passes.

$gene = new Bio::EnsEMBL::Gene;
print "ok 2\n";   
$tr   = new Bio::EnsEMBL::Transcript;
print "ok 3\n";   
$ex1   = new Bio::EnsEMBL::Exon;
$ex2   = new Bio::EnsEMBL::Exon;
print "ok 4\n";   


# to keep test quiet, initialise Gene/transcript at 20,20

#$gene->start(20);
#$gene->end(20);
#$tr->start(20);
#$tr->end(20);


$ex1->start(10);
$ex1->end(20);
$ex1->strand(1);

$ex2->start(40);
$ex2->end(50);
$ex2->strand(1);


$tr->add_Exon($ex1);
$tr->add_Exon($ex2);

$gene->add_Transcript($tr);

print "ok 5\n";

if( $gene->start != 10 || $gene->end != 50 ) {
    print "not ok 6\n";
} else {
    print "ok 6\n";
}

$count = 0;
foreach $tr ( $gene->each_Transcript() ) {
	foreach $x ( $tr->each_Exon() ) {
		$count++;
		}
	}

$x =0; # stop it whining...


if( $count != 2 ) {
    print "not ok 7\n";
} else {
    print "ok 7\n";
}

