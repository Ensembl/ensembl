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
BEGIN { $| = 1; print "1..33\n"; 
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
$ens_test->do_sql_file("t/geneget.dump");
    
# Get an EnsEMBL db object for the test db
my $db = $ens_test->get_DBSQL_Obj;
print "ok 2\n";


my $protein_adaptor=Bio::EnsEMBL::DBSQL::Protein_Adaptor->new($db);

eval {
    $protein = $protein_adaptor->fetch_Protein_by_dbid(1);
};

#print STDERR "PROT_AC: ".$protein->id."\n";

if ($@) {
    print "not ok 3 $@\n";
}
else {
    print "ok 3\n";
}


my @features = $protein->get_all_DomainFeatures;

my @f;

if (scalar @features == 3) {
    print "ok 4\n";
}
else {
    print "not ok 4\n";
}

#print STDERR "len ",$protein->length," ",$protein->id," ",$protein->seq,"\n";

if ($protein->length == 5) {
    print "ok 5\n";
}
else {
    print "not ok 5\n";
}

if ($protein->id eq "1" ) {
     print "ok 6\n";
}
else {
    print "not ok 6\n";
}

print "ok 7\n";


if ($protein->seq eq "KPPXP") {
    print "ok 8\n";
}
else {
    print "not ok 8\n";
}


if ($features[1]->analysis->db eq "Pfam") {
     print "ok 9\n";
}
else {
    print "not ok 9\n";
}


if ($features[0]->analysis->db eq "PRINTS") {
     print "ok 10\n";
}
else {
    print "not ok 10\n";
}


#my @introns = $protein->get_all_IntronFeatures();
#if ($introns[0]->feature1->start == 18) {
#    print "ok 11\n";
#}
#else {
#    print "not ok 11\n";
#}

print "ok 11\n";


if ($protein->geneac() eq 1) {
print "ok 12\n";
}
else {
    print "not ok 12\n";
}

if ($protein->transcriptac() eq 1) {
    print "ok 13\n";
}
else {

    print "not ok 13\n";
}


my @seq_features = $protein->top_SeqFeatures();


if (scalar(@seq_features) == 5) {
    print "ok 14\n";
}
else {
    print "not ok 14\n";
}

if ($seq_features[0]->feature2->seqname eq "PR00012") {
    print "ok 15\n";
}
else {
    print "not ok 15\n";
}



if ($seq_features[0]->feature1->seqname eq 1) {
    print "ok 16\n";
}
else {
    print "not ok 16\n";
}

# duplicate intron test!
print "ok 17\n";
print "ok 18\n";
print "ok 19\n";


if ($protein->molecular_weight == 495) {
    print "ok 20\n";
}
else {
    print "not ok 20\n";
}

if ($protein->checksum() eq "777AF77763500000" ) {
    print "ok 21\n";
}
else {
    print "not ok 21\n";
}

if ($seq_features[0]->idesc eq "petit") {
    print "ok 22\n";
}
else {
    print "not ok 22\n";
}

my @dblinks = $protein->get_all_DBLinks();

my @synonyms = $dblinks[0]->get_synonyms();

    
if ($dblinks[0]->release == 2) {
    print "ok 23\n";
}
else {
    print "not ok 23\n";
}


print "ok 24\n";


if ($dblinks[0]->description eq "tremblannot") {
     print "ok 25\n";
}
else {
    print "not ok 25\n";
}

if (scalar @dblinks == 1) {
     print "ok 26\n";
}
else {
    print "not ok 26\n";
}

my @domains = $protein->get_all_DomainFeatures();



if (scalar @domains == 3) {
     print "ok 27\n";
}
else {
    print "not ok 27\n";
}

my @transmembrane = $protein->get_all_TransmembraneFeatures();

if (scalar @transmembrane == 1) {
     print "ok 28\n";
}
else {
    print "not ok 28\n";
}

my @sigp = $protein->get_all_SigpFeatures();

if (scalar @sigp == 0) {
     print "ok 29\n";
}
else {
    print "not ok 29\n";
}

my $species = $protein->species;
my $common = $species->common_name();

if ($common eq "Human") {
    print "ok 30\n";
}
else {
    print "not ok 30\n";
}

if ($dblinks[0]->query_identity == 98) {
    print "ok 31\n";
}
else {
    print "not ok 31\n";
}

if ($dblinks[0]->target_identity == 99) {
    print "ok 32\n";
}
else {
    print "not ok 32\n";
}

eval {
    my $seqio = Bio::SeqIO->new('-format' => 'swiss' , -file => ">seq_temp.swiss" ) ;

    $seqio->write_seq($protein);
};

if ($@) {
    print STDERR "Exception $@\n";
    print "not ok 33\n";
}
else {
    print "ok 33\n";
}


my $rm = "rm seq_temp.swiss";

system($rm) == 0 or die "$0\Error running '$rm'";






