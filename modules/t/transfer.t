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

$conf{'donor'} = 'testdonor';
$conf{'recipient'} = 'testrecipient';
$conf{'mysqladmin'} = '/mysql/current/bin/mysqladmin';
$conf{'mysql'} = '/mysql/current/bin/mysql';
$conf{'user'}  = 'ensembl';
$conf{'update'} = '../scripts/update_list_chunk.pl';
$conf{'perl'} = 'perl';

if ( -e 't/transfer.conf' ) {
   print STDERR "Reading configuration from transfer.conf\n";
   open(C,"t/transfer.conf");
   while(<C>) {
       my ($key,$value) = split;
       $conf{$key} = $value;
   }
} else {
    print STDERR "Using default values\n";
    foreach $key ( keys %conf ) {
       print STDERR " $key $conf{$key}\n";
       }
    print STDERR "\nPlease use a file t/transfer.conf to alter these values if the test fails\nFile is written <key> <value> syntax\n\n";
}

$nuser = $conf{user};

my $create_donor = "$conf{mysqladmin} -u ".$nuser." create $conf{donor}";
my $create_recipient = "$conf{mysqladmin} -u ".$nuser." create $conf{recipient}";
system($create_donor) == 0 or die "$0\nError running '$create_donor' : $!";
system($create_recipient) == 0 or die "$0\nError running '$create_recipient' : $!";

print "ok 2\n";    #Databases created successfuly

#Initialising databases
my $init_donor = "$conf{mysql} -u ".$nuser." $conf{donor} < ../sql/table.sql";
my $init_recipient = "$conf{mysql} -u ".$nuser." $conf{recipient} < ../sql/table.sql";
system($init_donor) == 0 or die "$0\nError running '$init_donor' : $!";
system($init_recipient) == 0 or die "$0\nError running '$init_recipient' : $!";

print "ok 3\n";

#Suck test data into donor
print STDERR "Inserting test data in test donor db... this will take a while...\n";
my $suck_data = "$conf{mysql} -u ".$nuser." $conf{donor} < t/donor.dump";
system($suck_data) == 0 or die "$0\nError running '$suck_data' : $!";

print "ok 4\n";

#Insert values in meta table of recipient
my $meta= "echo \"insert into meta (donor_database_locator) values('Bio::EnsEMBL::DBSQL::Obj/host=localhost;port=410000;dbname=$conf{donor};user=$nuser;pass=');\" | $conf{mysql} -u $nuser $conf{recipient}";
system($meta) == 0 or die "$0\nError running '$meta' : $!";

print "ok 5\n";

#Update recipient from donor
print STDERR "Running an update from the donor to the recipient\n";
my $update="$conf{perl} $conf{update} -thost localhost -tdbname $conf{recipient} -tdbuser $nuser";
system($update) == 0 or die "$0\nError running '$update' : $!";

print "ok 6\n";

#Trying manual transfer of one gene from donor to recipient
my $dbtype = 'rdb';
my $host   = 'localhost';
my $port   = '410000';
my $dbname = $conf{'donor'};
my $dbuser = $nuser;
my $dbpass = undef;
my $module = 'Bio::EnsEMBL::DBSQL::Obj';

my $locator = "$module/host=$host;port=$port;dbname=$dbname;user=$dbuser;pass=$dbpass";
$don_db =  Bio::EnsEMBL::DBLoader->new($locator);

my $dbname = $conf{'recipient'};
my $locator = "$module/host=$host;port=$port;dbname=$dbname;user=$dbuser;pass=$dbpass";
$rec_db =  Bio::EnsEMBL::DBLoader->new($locator);

$gene=$don_db->get_Gene('ENSG00000019482');

$rec_db->delete_Gene($gene->id);
$rec_db->write_Gene($gene);
print "ok 7\n";

END {
    my $drop_donor = "echo \"y\" | $conf{mysqladmin} -u ".$nuser." drop $conf{donor}";
    my $drop_recipient = "echo \"y\" | $conf{mysqladmin} -u ".$nuser." drop $conf{recipient}";
    system($drop_donor) == 0 or die "$0\nError running '$drop_donor' : $!";
    system($drop_recipient) == 0 or die "$0\nError running '$drop_recipient' : $!";
}




