#!/usr/local/bin/perl

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB::CrossCloneMap;
use Bio::EnsEMBL::DBSQL::CrossMatchDBAdaptor;
use Bio::EnsEMBL::DBLoader;
use Getopt::Long;

$| = 1;


my $dbtype = 'rdb';
my $host   = 'ecs1a';
my $port   = '410000';
my $dbname = 'cross_oct07';
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
my $clone=shift @ARGV;

my $locator = "$module/host=$host;port=$port;dbname=$dbname;user=$dbuser;pass=$dbpass";

print STDERR "Using $locator for crossmatch db\n";
my $crossdb =  Bio::EnsEMBL::DBLoader->new($locator);
#my @clones=$crossdb->get_clonelist();

#foreach my $clone (@clones) {
print STDERR "Running mapping for clone $clone\n";
my $crossmap = Bio::EnsEMBL::Pipeline::RunnableDB::CrossCloneMap->new(-crossdb=>$crossdb, -score =>1000);
$crossmap->fetch_input($clone);

#$SIG{ALRM} = sub { die "timeout"};
#eval {
#    alarm(3600);
$crossmap->run;
#    alarm (0);
#};
#if ($@) {
#    if ($@ =~ /timeout/) {
#	die("EXCEPTION: Crossmatch for clone $clone timed out! Exiting");
#    }
#    else {
#	die("Died because of $@");
#    }
#}
print STDERR "Writing output for clone $clone\n";
$crossmap->write_output;



