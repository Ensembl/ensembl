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

use Bio::EnsEMBL::DB::Obj;
use Bio::EnsEMBL::DB::Clone;
use Bio::EnsEMBL::DB::Contig;

$loaded = 1;
print "ok 1\n";    # 1st test passes.

my $db;
eval {
    $db = new Bio::EnsEMBL::DB::Obj( -user => 'root', -db => 'pog' , -host => 'croc' );
};
if( $@ =~ /connect/ ) {
    print STDERR "Could not connect to DB. Skipping test\n";
    print "ok 2\n";
    print "ok 3\n";
    print "ok 4\n";
    print "ok 5\n";
    exit(0);
}


print "ok 2\n";
my $clone  = $db->get_Clone("dJ382I10");
print "ok 3\n";


my @contigs = $clone->get_all_Contigs();

my $contig = $db->get_Contig($contigs[0]->id);
print "ok 4\n";

foreach $sf ( $contig->get_all_SeqFeatures ) {
    # $sf is Bio::SeqFeature::Generic object.
    if( ! $sf->isa("Bio::SeqFeatureI") ) {
      print "not ok 5\n";
      exit(1);
    }
}

print "ok 5\n";

