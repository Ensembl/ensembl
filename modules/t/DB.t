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
BEGIN { $| = 1; print "1..11\n"; 
	use vars qw($loaded); }
END {print "not ok 1\n" unless $loaded;}

use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DBSQL::Clone;
use Bio::EnsEMBL::DBSQL::RawContig;
use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::DB::VirtualContig;

$loaded = 1;
print "ok 1\n";    # 1st test passes.

my $db;

sub skip_tests {
    print "ok 2\n";
    print "ok 3\n";
    print "ok 4\n";
    print "ok 5\n";
    exit(0);
}

open(FILE,"t/locator") || do {
		       print STDERR "Could not open locator string in t/locator\nYou need a database to test against in t/locator\nSee t/locator.example for an example locator\n";
		       print STDERR "\nDeliberately skipping tests\n";
		       &skip_tests();
		       exit(0);
		       };

$locator = <FILE>;
chomp $locator;		       
eval {
    $db = Bio::EnsEMBL::DBLoader->new($locator);	
};

if( $@  ) {
    print STDERR "Could not connect to database in locator [$locator]\nCompile went ok. Locator string looks incorrect\nDeliberately skipping test\n";
    &skip_tests();
}


print "ok 2\n";
@cloneids =  $db->get_all_Clone_id();
my $clone  = $db->get_Clone($cloneids[0]);

# check clone stuff.
$discard = $clone->htg_phase();
$discard = $clone->embl_id();
$discard = $clone->version();
$discard = $clone->embl_version();
print "ok 3\n";


my @contigs = $clone->get_all_Contigs();
my $contig = $db->get_Contig($contigs[0]->id);
print "ok 4\n";

@repeats = $contig->get_all_RepeatFeatures();
@repeats = ();
print "ok 5\n";

@simil   = $contig->get_all_SimilarityFeatures();
@simil = ();
print "ok 6\n";

foreach $gene ( $clone->get_all_Genes() ) {
    if( ! $gene->isa("Bio::EnsEMBL::Gene") ) {
      print "not ok 7\n";
      exit(1);
    }
}

print "ok 7\n";


eval {
    $contig = $db->get_Contig('test-contig-1');
};

if( $@ ) {
    print STDERR "Does not have contig test, cannot test overlap\n$@\n";
    print "ok 8\n";
    print "ok 9\n";
} else {
    $overlap = $contig->get_right_overlap();
    print "ok 8\n";


    $contig2 = $overlap->sister;
    $lefto = $contig2->get_left_overlap();
    
    print "ok 9\n";


    $vc = Bio::EnsEMBL::DB::VirtualContig->new( -focus => $contig2,
						-focusposition => 3,
						-ori => 1,
						-left => 6,
						-right => 9 );

    $vc->_dump_map();

    $seq = $vc->seq();
    if( $seq->isa('Bio::PrimarySeqI') ) {
	print "ok 10\n";
	print STDERR "Seq is ".$seq->seq."\n";
    } else {
	print "not ok 10\n";
    }

    @sf = $vc->get_all_SeqFeatures();
    foreach $sf ( @sf ) {
	print STDERR "Feature starts on ",$sf->start," ends ",$sf->end, "\n";
    }


    print "ok 11\n";

}
    





