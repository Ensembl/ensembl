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
BEGIN { $| = 1; print "1..7\n"; 
	use vars qw($loaded); }
END {print "not ok 1\n" unless $loaded;}

use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DBLoader;
$loaded=1;
print "ok \n";    # 1st test passed, loaded needed modules

#Creating test donor and recipient databases

print STDERR "Please insert read/write user (not password protected, e.g. ensembl): ";
$nuser=<STDIN>;
chop $nuser;

my $create_donor = "mysqladmin -u ".$nuser." create donor";
my $create_recipient = "mysqladmin -u ".$nuser." create recipient";
system($create_donor) == 0 or die "$0\nError running '$create_donor' : $!";
system($create_recipient) == 0 or die "$0\nError running '$create_recipient' : $!";

print "ok 2\n";    #Databases created successfuly

#Initialising databases
my $init_donor = "mysql -u ".$nuser." donor < ../sql/table.sql";
my $init_recipient = "mysql -u ".$nuser." recipient < ../sql/table.sql";
system($init_donor) == 0 or die "$0\nError running '$init_donor' : $!";
system($init_recipient) == 0 or die "$0\nError running '$init_recipient' : $!";

print "ok 3\n";

#Suck test data into donor
print STDERR "Inserting test data in test donor db... this will take a while...\n";
my $suck_data = "mysql -u ".$nuser." donor < donor.dump";
system($suck_data) == 0 or die "$0\nError running '$suck_data' : $!";

print "ok 4\n";

#Insert values in meta table of recipient
my $meta= "echo \"insert into meta (donor_database_locator) values('Bio::EnsEMBL::DBSQL::Obj/host=localhost;port=410000;dbname=donor;user=ensembl;pass=');\" | mysql -u $nuser recipient";
system($meta) == 0 or die "$0\nError running '$meta' : $!";

print "ok 5\n";

#Update recipient from donor
my $update="perl ../../scripts/update_list_chunk.pl -thost localhost -tdbname recipient -tdbuser $nuser";
#system($update) == 0 or die "$0\nError running '$update' : $!";

print "ok 6\n";

my $drop_donor = "echo \"y\" | mysqladmin -u ".$nuser." drop donor";
my $drop_recipient = "echo \"y\" | mysqladmin -u ".$nuser." drop recipient";
system($drop_donor) == 0 or die "$0\nError running '$drop_donor' : $!";
system($drop_recipient) == 0 or die "$0\nError running '$drop_recipient' : $!";
print "ok 7\n";



