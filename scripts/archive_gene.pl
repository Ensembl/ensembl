#!usr/local/bin/perl

use strict;

use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DBLoader;
use Getopt::Long;
use Bio::Seq;

my $tdbtype = 'rdb';
my $thost   = 'obi-wan.sanger.ac.uk';
my $tport   = '410000';
my $tdbname = 'ensembl';
my $tdbuser = 'root';
my $tpass = undef;
my $adbname = 'ens_archive';
my $module = "Bio::EnsEMBL::DBSQL::Obj";

my $dbtype = 'rdb';
my $host   = 'localhost';
my $port   = '410000';
my $dbname = 'archive';
my $dbuser = 'root';
my $arcmodule = 'Bio::EnsEMBL::DBArchive::Obj';
my $help;

&GetOptions( 
	     'tdbtype:s'  => \$tdbtype,
	     'thost:s'    => \$thost,
	     'tport:n'    => \$tport,
	     'tdbname:s'  => \$tdbname,
	     'tdbuser=s'  => \$tdbuser,
	     'tpass:s'    => \$tpass,
	     'dbtype:s'   => \$dbtype,
	     'host:s'     => \$host,
	     'port:n'     => \$port,
	     'dbname:s'   => \$dbname,
	     'dbuser:s'   => \$dbuser,
	     'module:s'   => \$module,
	     'h|help'     => \$help
	     );

my @genes = @ARGV;

my $arclocator = "$arcmodule/host=$host;port=$port;dbname=$dbname;user=$dbuser;pass=$tpass;debug=10";
my $arcdb =  Bio::EnsEMBL::DBLoader->new($arclocator);

my $to_locator = "$module/host=$thost;port=$tport;dbname=$tdbname;user=$tdbuser;pass=$tpass";
my $tdb = new Bio::EnsEMBL::DBLoader($to_locator);

foreach my $geneid (@genes){
    print STDERR "Archiving $geneid\n";
    my $gene=$tdb->get_Gene($geneid);
    $tdb->archive_Gene($gene,$arcdb);
}


