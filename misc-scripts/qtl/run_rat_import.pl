
#
# Program to download latest QTLs from RGD, delete old qtls in the database and
# import the new ones.
#

use strict;

use Getopt::Long;
use DBI;



my $RGD_URL    = 'ftp://rgd.mcw.edu/pub/data_release/QTLS';

my ($user, $pass, $port, $dbname, $host);

GetOptions( "host=s", \$host,
	    "user=s", \$user,
	    "pass=s", \$pass,
	    "port=i", \$port,
	    "dbname=s", \$dbname,);

$port ||= 3306;

if( !$host || !$dbname ) {
  usage();
}

my $dsn = "DBI:mysql:host=$host;dbname=$dbname;port=$port";
my $db = DBI->connect( $dsn, $user, $pass, {'RaiseError' => 1})
  or die("Could not connect to database $dbname $!");

print STDERR "Downloading RGD data\n";
system("/usr/local/bin/wget $RGD_URL -O /tmp/rgd.txt");

die("Download failed.") if(!-e '/tmp/rgd.txt');

print STDERR "Deleting existing QTL data\n";

$db->do('DELETE FROM qtl');
$db->do('DELETE FROM qtl_feature');
$db->do('DELETE FROM qtl_synonym');

print STDERR "Importing QTLs\n";

system("/usr/local/ensembl/bin/perl rat_qtl_import.pl -host $host " .
       " -user $user -pass $pass -port $port -dbname $dbname -verbose " .
       " -qtlfile /tmp/rgd.txt > /tmp/qtl.sql");

unlink('/tmp/rgd.txt');

system("cat /tmp/qtl.sql | /usr/local/bin/mysql -h $host -u $user -p$pass " .
       " -P$port $dbname");

unlink('/tmp/qtl.sql');

print STDERR "Generating QTL Features\n";

system("/usr/local/ensembl/bin/perl qtl_feature_calculation.pl " .
       "-host $host -user $user -pass $pass -port $port -dbname $dbname " .
       "-verbose > /tmp/qtl_feature.txt");

system("/usr/local/bin/mysqlimport -h $host -u $user -p$pass -P$port " .
       " $dbname /tmp/qtl_feature.txt");

unlink('/tmp/qtl_feature.txt');


$db->disconnect();

print STDERR "All done\n";


sub usage {
  print STDERR "usage:\n" .
               "  perl run_rat_import.pl -host <host> -user <user> " .
               "-dbname <dbname> [-pass <pass>] [-port <port>] \n" .
  exit;
}
