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
BEGIN { $| = 1; print "1..29\n"; 
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

if (scalar @features == 3) {
    print "ok 5\n";
}
else {
    print "not ok 5\n";
}

if ($protein->length == 46) {
    print "ok 6\n";
}
else {
    print "not ok 6\n";
}

if ($protein->id eq "ENSP00000216167") {
     print "ok 7\n";
}
else {
    print "not ok 7\n";
}

my @dates = $protein->get_dates();

if (scalar @dates == 2) {
     print "ok 8\n";
}
else {
    print "not ok 8\n";
}

my @dblinks = $protein->annotation->each_DBLink();

my @synonyms = $dblinks[0]->get_synonyms();

if ($dblinks[0]->url_pattern eq "http//") {
    print "ok 9\n";
}
else {
    print "not ok 9\n";
}

    
if ($dblinks[0]->release == 23) {
    print "ok 10\n";
}
else {
    print "not ok 10\n";
}

if ($synonyms[0] eq "EMBL_SYNONYME1") {
    print "ok 11\n";
}
else {
    print "not ok 11\n";
}

if ($dblinks[1]->description eq "tremblannot") {
     print "ok 12\n";
}
else {
    print "not ok 12\n";
}

if (scalar @dblinks == 7) {
     print "ok 13\n";
}
else {
    print "not ok 13\n";
}

if ($protein->seq eq "RNSKRTLCMNNLFPHYRQKNPRLLREPSDFLHLKSVKSSCFLLPYP") {
    print "ok 14\n";
}
else {
    print "not ok 14\n";
}

if ($features[0]->analysis->db eq "Pfam") {
     print "ok 15\n";
}
else {
    print "not ok 15\n";
}



if ($features[1]->analysis->db eq "PRINTS") {
     print "ok 16\n";
}
else {
    print "not ok 16\n";
}

my @introns = $protein->each_Intron_feature();

if ($introns[0]->feature1->start == 18) {
    print "ok 17\n";
}
else {
    print "not ok 17\n";
}

my @features2 = $protein->each_Protein_feature();

if (scalar @features2 == 3) {
    print "ok 18\n";
}
else {
    print "not ok 18\n";
}

if ($protein->geneac() eq "ENSG00000100331") {
print "ok 19\n";
}
else {
    print "not ok 19\n";
}

if ($protein->transcriptac() eq "ENST00000216167") {
    print "ok 20\n";
}
else {

    print "not ok 20\n";
}

my @seq_features = $protein->top_SeqFeatures();



if (scalar(@seq_features) == 3) {
    print "ok 21\n";
}
else {
    print "not ok 21\n";
}

if ($seq_features[0]->feature2->seqname eq "PF00098") {
    print "ok 22\n";
}
else {
    print "not ok 22\n";
}



if ($seq_features[0]->feature1->seqname eq "ENSP00000216167") {
    print "ok 23\n";
}
else {
    print "not ok 23\n";
}

my @introns = $protein->each_Intron_feature();

if (scalar(@introns) == 1) {
    
    print "ok 24\n";
}
else {
    print "not ok 24\n";
}
    

if ($introns[0]->intron_position == 18) {
    
    print "ok 25\n";
}
else {
    print "not ok 25\n";
}


if ($introns[0]->intron_length == 94) {
    
    print "ok 26\n";
}
else {
    print "not ok 26\n";
}


if ($protein->molecular_weight == 5547) {
    print "ok 27\n";
}
else {
    print "not ok 27\n";
}

if ($protein->checksum() eq "D5C2BBEAA77A0FB6") {
    print "ok 28\n";
}
else {
    print "not ok 28\n";
}

if ($seq_features[0]->idesc eq "Zn-finger") {
    print "ok 29\n";
}
else {
    print "not ok 29\n";
}


my $rm = "rm seq_temp.swiss";

system($rm) == 0 or die "$0\Error running '$rm'";










