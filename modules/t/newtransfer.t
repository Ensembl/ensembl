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
BEGIN { $| = 1; print "1..3\n"; 
	use vars qw($loaded); }

END {print "not ok 1\n" unless $loaded;
     system "rm t/recipient.dump";}

use lib 't';
use EnsTestDB;
$loaded = 1;
print "ok 1\n";    # 1st test passes.
    
my $ens_donor = EnsTestDB->new();
    
# Load some data into the db
$ens_donor->do_sql_file("t/donor.dump");
    
# Get an EnsEMBL db object for the test db
my $fromdb = $ens_donor->get_DBSQL_Obj;
print "ok 2\n";    

my $ens_recipient = EnsTestDB->new();
open (TEMP,">t/recipient.dump");
my $meta= "insert into meta (donor_database_locator) values('Bio::EnsEMBL::DBSQL::Obj/host=".$ens_donor->host.";port=".$ens_donor->port.";dbname=".$ens_donor->dbname.";user=".$ens_donor->user.";pass=".$ens_donor->password."')";
print TEMP $meta;
$ens_recipient->do_sql_file("t/recipient.dump");
my $todb = $ens_recipient->get_DBSQL_Obj;
print "ok 3\n";

my @clones = $fromdb->get_Update_Obj->get_updated_Clone_id(1,time());
print "ok 4\n";


foreach my $clone_id (@clones) {
    
    my $don_clone=$fromdb->get_Clone($clone_id);
    $todb->write_Clone($don_clone);
    foreach my $gene ($don_clone->get_all_Genes('evidence')) {
	$todb->gene_Obj->write($gene);
    }
}
print "ok 5\n";











































