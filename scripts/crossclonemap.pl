#!/usr/local/bin/perl

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB::CrossCloneMap;
use Bio::EnsEMBL::DBSQL::CrossMatchDBAdaptor;
use Bio::EnsEMBL::DBLoader;
use Getopt::Long;

$| = 1;


my $dbtype = 'rdb';
my $host   = 'ecs1c';
my $port   = '410000';
my $dbname = 'cross110_bis';
my $dbuser = 'ensadmin';
my $dbpass = undef;
my $module = 'Bio::EnsEMBL::DBSQL::CrossMatchDBAdaptor';

&GetOptions ( 
	       'dbtype:s' => \$dbtype,
	       'host:s'   => \$host,
	       'port:n'   => \$port,
	       'dbname:s' => \$dbname, 
	       'dbuser:s' => \$dbuser,
	       'dbpass:s' => \$dbpass,
	       'module:s' => \$module,
	       );
my $locator = "$module/host=$host;port=$port;dbname=$dbname;user=$dbuser;pass=$dbpass";

print STDERR "Using $locator for crossmatch db\n";
my $crossdb =  Bio::EnsEMBL::DBLoader->new($locator);
my @clones=$crossdb->get_clonelist();

foreach my $clone (@clones) {
    print STDERR "Sending crossclonemap job for clone $clone to LSF queue\n";
    my $command = "bsub -o $clone.out -e $clone.err -E /work2/elia/src/scripts/echeck.pl /work2/elia/src/ensembl/scripts/clonemap.pl $clone";
    print STDERR "Command: $command\n";
    system($command);
}
$crossdb->DESTROY;


