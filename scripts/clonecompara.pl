#!/usr/local/bin/perl

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB::CrossComparer;
use Bio::EnsEMBL::DBLoader;
use Getopt::Long;

$| = 1;


my $dbtype = 'rdb';
my $host   = 'ensrv3';
my $port   = '';
my $dbname = 'homo_sapiens_core_110';
my $dbuser = 'ensro';
my $dbpass = '';
my $module = 'Bio::EnsEMBL::DBSQL::DBAdaptor';

my $dbtype2 = 'rdb';
my $host2   = 'ensrv3';
my $port2   = '';
my $dbname2 = 'homo_sapiens_core_110';
my $dbuser2 = 'ensro';
my $dbpass2 = '';
my $module2 = 'Bio::EnsEMBL::DBSQL::DBAdaptor';
my $contig;
my $contig2;

&GetOptions ( 
	      'dbtype:s' => \$dbtype,
	      'host:s'   => \$host,
	      'port:n'   => \$port,
	      'dbname:s' => \$dbname, 
	      'dbuser:s' => \$dbuser,
	      'dbpass:s' => \$dbpass,
	      'module:s' => \$module,
	      'dbtype2:s' => \$dbtype2,
	      'host2:s'   => \$host2,
	      'port2:n'   => \$port2,
	      'dbname2:s' => \$dbname2, 
	      'dbuser2:s' => \$dbuser2,
	      'dbpass2:s' => \$dbpass2,
	      'module2:s' => \$module2,
	      'contig1:s' => \$contig,
	      'contig2:s' => \$contig2
	      );
$contig || die("Need to specify contigs to compare (got $contig and $contig2)");
$contig2 || die("Need to specify contigs to compare (got $contig and $contig2)");

my $locator = "$module/host=$host;port=$port;dbname=$dbname;user=$dbuser;pass=$dbpass";
print STDERR "Using $locator for db 1\n";

my $locator2 = "$module2/host=$host2;port=$port2;dbname=$dbname2;user=$dbuser2;pass=$dbpass2";
print STDERR "Using $locator for db 2\n";

my $querydb =  Bio::EnsEMBL::DBLoader->new($locator);
my $targetdb =  Bio::EnsEMBL::DBLoader->new($locator2);

my $input_id = "$locator-$contig~$locator2-$contig2";

my $crosscomparer = Bio::EnsEMBL::Pipeline::RunnableDB::CrossComparer->new(
									 -dbobj => $querydb,
									 -input_id => $input_id,
									 -min_score => 100
									 );

$crosscomparer->fetch_input;
$crosscomparer->run;
$crosscomparer->write_output;




