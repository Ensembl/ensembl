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

#
#NB: Are not tested the following methods: get_Protein_annseq, write_all_Protein_features, write_Protein_feature.
#

## We start with some black magic to print on failure.
BEGIN { $| = 1; print "1..5\n"; 
	use vars qw($loaded); }
END {print "not ok 1\n" unless $loaded;}


use lib 't';
use EnsTestDB;
use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DBSQL::Protein_Adaptor;
use Bio::EnsEMBL::Protein;

$loaded = 1;
print "ok 1\n";    # 1st test passes.
    
my $ens_test = EnsTestDB->new();
    
# Load some data into the db
$ens_test->do_sql_file("t/protein.dump");
    
# Get an EnsEMBL db object for the test db
my $db = $ens_test->get_DBSQL_Obj;
print "ok 2\n";


 my $protein_adaptor=Bio::EnsEMBL::DBSQL::Protein_Adaptor->new($db);
eval {
    $protein = $protein_adaptor->fetch_Protein_by_dbid('ENSP00000216167');
};
if ($@) {
    print "not ok 3\n";
}
else {
    print "ok 3\n";
}

eval {
    my $seqio = Bio::SeqIO->new('-format' => 'swiss' , -file => ">seq_temp.swiss" ) ;

    $seqio->write_seq($protein);
};

if ($@) {
    print STDERR "Exception $@\n";
    print "not ok 4\n";
}
else {
    print "ok 4\n";
}


my @features = $protein->each_Protein_feature;

if (scalar @features == 2) {
    print "ok 5\n";
}
else {
    print "not ok 5\n";
}

my $rm = "rm seq_temp.swiss";

system($rm) == 0 or die "$0\Error running '$rm'";
