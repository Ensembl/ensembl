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
# This test does not test of all the components of DBEntryAdaptor, most of them are already tested in protein.t

## We start with some black magic to print on failure.
BEGIN { $| = 1; print "1..15\n"; 
	use vars qw($loaded); }
END {print "not ok 1\n" unless $loaded;}


use lib 't';
use EnsTestDB;
use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DBSQL::DBEntryAdaptor;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::IdentityXref;

$loaded = 1;
print "ok 1\n";    # 1st test passes.
    
my $ens_test = EnsTestDB->new();
    
# Load some data into the db
$ens_test->do_sql_file("t/DBEntryAdaptor.dump");
    
# Get an EnsEMBL db object for the test db
my $db = $ens_test->get_DBSQL_Obj;
print "ok 2\n";


my $dbentry_adaptor = Bio::EnsEMBL::DBSQL::DBEntryAdaptor->new($db);

my $dbentry = Bio::EnsEMBL::DBEntry->new 
    ( -primary_id => 123,
      -display_id => 'Ensembl_link',
      -version => 2,
      -release => 2344,
      -dbname =>  'ENSDBNAME');

eval {
    $dbentry_adaptor->store($dbentry, 10, 'translation');
};

if ($@) {
    print "not ok 3\n";
}

else {
    print "ok 3\n";
}

my $entry1;

eval {
    ($entry1) = $dbentry_adaptor->fetch_by_translation(10);
};

if ($@) {
    print "not ok 4\n";
}

else {
    print "ok 4\n";
}

if ($entry1->dbname eq "ENSDBNAME") {
    print "ok 5\n";
}

else {
    print "not ok 5\n";
}

if ($entry1->version == 2) {
    print "ok 6\n";
}

else {
    print "not ok 6\n";
}

my $dbentry2 = Bio::EnsEMBL::IdentityXref->new 
    ( -primary_id => 245,
      -display_id => 'Ensembl_link2',
      -version => 2,
      -release => 2344,
      -dbname =>  'ENSDBNAME_idt');

$dbentry2->target_identity(67);
$dbentry2->query_identity(100);

eval {
    $dbentry_adaptor->store($dbentry2, 11, 'translation');
};

if ($@) {
    print "not ok 7\n";
}

else {
    print "ok 7\n";
}

my $dbentry4;

eval {
    ($entry4) = $dbentry_adaptor->fetch_by_translation(11);
};

if ($@) {
    print "not ok 8\n";
}

else {
    print "ok 8\n";
}

if ($entry4->query_identity == 100) {
     print "ok 9\n";
}

else {
    print "not ok 9\n";
}

if ($entry4->target_identity == 67) {
     print "ok 10\n";
}

else {
    print "not ok 10\n";
}

my @matches = $dbentry_adaptor->fetch_by_gene(1);

if (scalar (@matches) == 1) {
    print "ok 11\n";
}

else {
    print "not ok 11\n";
}

if ($matches[0]->primary_id eq "Q9Ngene") {
    print "ok 12\n";
}

else {
    print "not ok 12\n";
}

my @matches2 = $dbentry_adaptor->fetch_by_translation(1);


if (scalar (@matches2) == 2) {
    print "ok 13\n";
}

else {
    print "not ok 13\n";
}

if ($matches2[0]->primary_id eq "AL365511") {
     print "ok 14\n";
}

else {
    print "not ok 14\n";
}

if ($matches2[1]->primary_id eq "Q9NPS7") {
     print "ok 15\n";
}

else {
    print "not ok 15\n";
}

