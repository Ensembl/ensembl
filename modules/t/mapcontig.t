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
BEGIN { $| = 1; print "1..9\n"; 
	use vars qw($loaded); }
END {print "not ok 1\n" unless $loaded;}

use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::Virtual::MapContig;
use Bio::EnsEMBL::Virtual::Map;

use lib 't';
use EnsTestDB;
$loaded = 1;
print "ok 1\n";    # 1st test passes.
    
my $ens_test = EnsTestDB->new();
    
# Load some data into the db
$ens_test->do_sql_file("t/mapcontig.dump");
    
# Get an EnsEMBL db object for the test db
my $db = $ens_test->get_DBSQL_Obj;
print "ok 2\n";    

my $contig = $db->get_Contig('contig1');

die "$0\nError fetching contig1 : $!" unless defined ($contig);

print "ok 3\n";


$mc = Bio::EnsEMBL::Virtual::MapContig->new(
	-rawcontig => $contig,
	-start => 400,
	-end   => 420,
	-rawcontig_start => 3,
	-orientation => 1
	);

print "ok 4\n";

if( $mc->start != 400 ||
    $mc->end   != 420 ||
    $mc->rawcontig_start != 3 ||
    $mc->rawcontig_end   != 23 ||
    $mc->orientation     != 1 ) {
    print STDERR "Map contig managed to get mangled\n";
    print "not ok 5\n";
} else {
    print "ok 5\n";
}


$map = Bio::EnsEMBL::Virtual::Map->new();

$map->create_MapContig($contig,400,420,3,1);

$mc = $map->get_MapContig_by_id($contig->id);


if( $mc->start != 400 ||
    $mc->end   != 420 ||
    $mc->rawcontig_start != 3 ||
    $mc->rawcontig_end   != 23 ||
    $mc->orientation     != 1 ) {
    print STDERR "Map contig managed to get mangled\n";
    print "not ok 6\n";
} else {
    print "ok 6\n";
}

# raw contig seq  is:
# AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT
$str = $mc->_actual_sequence_as_string;
$shouldbe='AAAAAAAACCCCCCCCCCGGG';
# if( $str ne 'AAAAAAAATTTTTTTTTAAAA' ) {
  if( $str ne $shouldbe ) {
    print "not ok 7\n";
    print STDERR "Seq $str, should be $shouldbe\n";
} else {
    print "ok 7\n";
}

# reverse strand:
$mc = Bio::EnsEMBL::Virtual::MapContig->new(
	-rawcontig => $contig,
	-start => 1001,
	-end   => 1010,
 	-rawcontig_start => 8,
	-orientation => -1,
	);
# raw contig seq  is:
# AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT
$str=$mc->_actual_sequence_as_string ;
$shouldbe='GGGGGGGTTT';
if ($str eq $shouldbe )  {
  print "ok 8\n";
} else {
  print "not ok 8\n";
  warn "Seq $str, should be $shouldbe\n";
}

# check missing args:
eval { 
  $mc = Bio::EnsEMBL::Virtual::MapContig->new(
	-rawcontig => $contig,
	-start => 230,
 	-rawcontig_start => 3,
	-orientation => 1
	);
}; 
if ($@) {
   print "ok 8\n";
} else {
   print "not ok 8\n";
   warn "expected exception on missing arguments";
}

eval {
$mc = Bio::EnsEMBL::Virtual::MapContig->new(
	-rawcontig => $contig,
	-start => 30,
	-end   => 10,
 	-rawcontig_start => 3,
	-orientation => 1
	);
};
if ( $@ ) {
   print "ok 9\n";
} else {
   print "not ok 9\n";
   warn "expected exception on start > end ";
}
