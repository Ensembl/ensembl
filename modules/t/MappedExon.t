## Bioperl Test Harness Script for Modules
##


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

use Bio::EnsEMBL::MappedExon;

$loaded = 1;
print "ok 1\n";    # 1st test passes.


my $exon = Bio::EnsEMBL::MappedExon->new();
$exon->start(10);
$exon->end(20);
$exon->strand(1);
$exon->contig_id('AC00013.1');

if (!defined($exon)) {
    print "not ok 2\n";
} else {
    print "ok 2\n";
} 

$exon->has_identical_sequence(0);

print	"ok 3\n";

$exon->has_identical_sequence;

print	"ok 4\n";

