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

use Bio::EnsEMBL::PerlDB::Clone;
use Bio::EnsEMBL::PerlDB::Contig;
use Bio::EnsEMBL::PerlDB::Obj;

use Bio::Seq;
use Bio::EnsEMBL::Gene;

$loaded = 1;
print "ok 1\n";    # 1st test passes.


# make some new objects

$seq = Bio::Seq->new(-seq => 'AAATAATATATAATATATATAT' , -id => 'dummy');

$gene = new Bio::EnsEMBL::Gene;
$tr   = new Bio::EnsEMBL::Transcript;
$ex1   = new Bio::EnsEMBL::Exon;
$ex2   = new Bio::EnsEMBL::Exon;

$ex1->start(10);
$ex1->end(20);
$ex1->strand(1);

$ex2->start(40);
$ex2->end(50);
$ex2->strand(1);


$tr->add_Exon($ex1);
$tr->add_Exon($ex2);

$gene->add_Transcript($tr);

# make a contig object

$contig = Bio::EnsEMBL::PerlDB::Contig->new();

# add sequence and genes

$contig->id("contig_id");
$contig->seq($seq);
$contig->add_Gene($gene);

print "ok 2\n";

# make a clone object

$clone = Bio::EnsEMBL::PerlDB::Clone->new();

# add contig object

$clone->add_Contig($contig);

print "ok 3\n";

# get stuff out.

foreach $contig ( $clone->get_all_Contigs ) {
	if( $contig->seq->id() ne 'dummy' ) {
	    print "not ok 4\n";
	} else {
	    print "ok 4\n";
        }

}

# make a Obj for db.

$db = Bio::EnsEMBL::PerlDB::Obj->new();

$db->write_Contig($contig);

$c = $db->get_Contig("contig_id");
$c =0;
print "ok 5\n";
