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
BEGIN { $| = 1; print "1..2\n"; 
	use vars qw($loaded); }

END {print "not ok 1\n" unless $loaded;}

use Bio::EnsEMBL::Ghost;
print "ok 1\n";    # 1st test passes.

$ghost = new Bio::EnsEMBL::Ghost;
$ghost->id('bollocks');
$ghost->version(2);
$ghost->obj_type('Gene');
$ghost->deleted(924455545);
$ghost->add_pointer_id('this');
$ghost->add_pointer_id('that');
@ids = $ghost->each_pointer_id();
if( scalar @ids != 2 ) {
  print "not ok 2\n"; 
} else {
  print "ok 2\n";
}





	  



