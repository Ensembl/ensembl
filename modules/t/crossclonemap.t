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

END { print "not ok 1\n" unless $loaded; }


use Bio::EnsEMBL::Pipeline::RunnableDB::CrossCloneMap;
use Bio::EnsEMBL::DBSQL::CrossMatchDBAdaptor;
#use Bio::EnsEMBL::DBLoader;
use lib 't';
use EnsTestDB;
print "ok 1\n";

#Create the old and new comparison databases;
my $old_testdb = EnsTestDB->new();
$old_testdb->do_sql_file("t/old_cross.dump");
my $old_db = $old_testdb->get_DBSQL_Obj;
print "ok 2\n";    

my $new_testdb = EnsTestDB->new();
$new_testdb->do_sql_file("t/new_cross.dump");
my $new_db = $new_testdb->get_DBSQL_Obj;
print "ok 3\n"; 

#Create the crossmatch database
my $test_db = EnsTestDB->new();
# Get an EnsEMBL db object for the test db
open (TEMP,">t/crossmatch.dump");
my $crloc = "INSERT INTO dblocation (olddatabase,newdatabase) VALUES ('Bio::EnsEMBL::DBSQL::Obj/host=".$old_testdb->host.";port=".$old_testdb->port.";dbname=".$old_testdb->dbname.";user=".$old_testdb->user.";pass=".$old_testdb->password."','Bio::EnsEMBL::DBSQL::Obj/host=".$new_testdb->host.";port=".$new_testdb->port.";dbname=".$new_testdb->dbname.";user=".$new_testdb->user.";pass=".$new_testdb->password."');\n";
print TEMP $crloc;
close (TEMP);

$test_db->do_sql_file('../sql/crossmatch.sql','t/crossmatch.dump');
my $host = $test_db->host;
my $dbname = $test_db->dbname;
my $user = $test_db->user;
my $crossdb = Bio::EnsEMBL::DBSQL::CrossMatchDBAdaptor->new( -host => $host, -dbname => $dbname, -user => $user );
print "ok 4\n"; 
system ('rm t/crossmatch.dump');

my $crossmap = Bio::EnsEMBL::Pipeline::RunnableDB::CrossCloneMap->new($crossdb);
print "ok 5\n";
$crossmap->fetch_input('crosstest_1');
print "ok 6\n";
$crossmap->run;
print "ok 7\n";

