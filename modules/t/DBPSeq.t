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
BEGIN { $| = 1; print "1..9\n"; 
	use vars qw($loaded); }
END {print "not ok 1\n" unless $loaded;}

use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DBLoader;
$loaded=1;
print "ok \n";    # 1st test passed, loaded needed modules

$conf{'mysqladmin'} = '/mysql/current/bin/mysqladmin';
$conf{'mysql'} = '/mysql/current/bin/mysql';
$conf{'user'}  = 'root';
$conf{'database'} = 'ensembl07';
$conf{'contig'} = 'AP000869.00011';

if ( -e 't/DBPSeq.conf' ) {
   print STDERR "Reading configuration from DBPSeq.conf\n";
   open(C,"t/DBPSeq.conf");
   while(<C>) {
       my ($key,$value) = split;
       $conf{$key} = $value;
   }
} else {
    print STDERR "Using default values\n";
    foreach $key ( keys %conf ) {
	print STDERR " $key $conf{$key}\n";
    }
    print STDERR "\nPlease use a file t/transf_to_vc.conf to alter these values if the test fails\nFile is written <key> <value> syntax\n\n";
}

my $dbtype = 'rdb';
my $host   = 'localhost';
my $port   = '410000';
my $dbname = $conf{'database'};
my $dbuser = $conf{'user'};
my $dbpass = undef;
my $module = 'Bio::EnsEMBL::DBSQL::Obj';

#Connect to local ensembl db
my $locator = "$module/host=$host;port=$port;dbname=$dbname;user=$dbuser;pass=$dbpass";
my $db =  Bio::EnsEMBL::DBLoader->new($locator);
print "ok 2\n";

#Get Contig
my $contig=$db->get_Contig($conf{'contig'});
print "ok 3\n";

#Get DBPrimarySeq for this contig
my $dbseq = $contig->DB_primary_seq();
print "ok 4\n";

#Check that DBPrimarySeq methods work
my $contig_id = $dbseq->contig_id;
print "ok 5\n";
#print STDERR "Internal id of $conf{'contig'} is $contig_id\n";

my $dna_id = $dbseq->dna_id;
print "ok 6\n";
#print STDERR "Dna id of $conf{'contig'} is $dna_id\n";

my $seq = $dbseq->seq();
print "ok 7\n";
#print STDERR "Seq of $conf{'contig'} is $seq\n\n";

my $subseq = $dbseq->subseq(100,200);
print "ok 8\n";
#print STDERR "Subseq from 100 to 200 is $subseq\n\n";

my $length = $dbseq->length();
print "ok 9\n";
#print STDERR "Length of dna of $conf{'contig'} is $length\n";


				      
