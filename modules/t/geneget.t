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
BEGIN { $| = 1; print "1..8\n"; 
	use vars qw($loaded); }
END {print "not ok 1\n" unless $loaded;}

use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DBLoader;
$loaded=1;
print "ok \n";    # 1st test passed, loaded needed modules

#Creating test overlap database

$conf{'overlap'}      = 'testoverlap';
$conf{'mysqladmin'} = '/mysql/current/bin/mysqladmin';
$conf{'mysql'}      = '/mysql/current/bin/mysql';
$conf{'user'}       = 'root';
$conf{'perl'}       = 'perl';

if ( -e 't/mysql.conf' ) {
  print STDERR "Reading configuration from overlap.conf\n";
  open(C,"t/mysql.conf");
  while(<C>) {
    my ($key,$value) = split;
    $conf{$key} = $value;
  }
} else {
  print STDERR "Using default values\n";
  foreach $key ( keys %conf ) {
    print STDERR " $key $conf{$key}\n";
  }
  print STDERR "\nPlease use a file t/mysql.conf to alter these values if the test fails\nFile is written <key> <value> syntax\n\n";
}

$nuser = $conf{user};

my $create_overlap        = "$conf{mysqladmin} -u ".$nuser." create $conf{overlap}";

system($create_overlap)   == 0 or die "$0\nError running '$create_overlap' : $!";

print "ok 2\n";    #Databases created successfuly

#Initialising databases
my $init_overlap        = "$conf{mysql} -u ".$nuser." $conf{overlap} < ../sql/table.sql";

system($init_overlap)     == 0 or die "$0\nError running '$init_overlap' : $!";

print "ok 3\n";

#Suck test data into db
print STDERR "Inserting test data in test overlap db...\n";
my $suck_data      = "$conf{mysql} -u ".$nuser." $conf{overlap} < t/geneget.dump";
system($suck_data) == 0 or die "$0\nError running '$suck_data' : $!";

print "ok 4\n";

# Connect to test db

my $db             = new Bio::EnsEMBL::DBSQL::Obj(-host   => 'localhost',
						  -user   => $conf{user},
						  -dbname => $conf{overlap}
						 );

die "$0\nError connecting to database : $!" unless defined($db);

print "ok 5\n";

$gene = $db->get_Gene('gene-id-1');

print "ok 6\n";

@dblink = $gene->each_DBLink();
$dbl = shift @dblink;

if( !defined $dbl || $dbl->database ne 'swissprot' ) {
    print "not ok 7\n";
} else {
  print "ok 7\n";
}

@trans = $gene->each_Transcript();
$trans = shift @trans;


@dblink = $trans->each_DBLink();
$dbl = shift @dblink;

if( !defined $dbl || $dbl->database ne 'embl' ) {
    print "not ok 8\n";
}    else {
  print "ok 8\n";
}


END {
    my $drop_overlap        = "echo \"y\" | $conf{mysqladmin} -u ".$nuser." drop $conf{overlap}";
   system($drop_overlap)     == 0 or die "$0\nError running '$drop_overlap' : $!";
}




