#!/usr/local/bin/perl

use strict;

use Getopt::Long;
use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DBLoader;

# global defaults
my $host = 'localhost';
my $dbuser = 'root';
my $dbname = 'newschema';
my $dbpass = undef;
my $module = 'Bio::EnsEMBL::DBSQL::Obj';
my $port = '410000';
my $usefile = 0;

&GetOptions( 
	     'dbuser:s'  => \$dbuser,
	     'dbpass:s'  => \$dbpass,
	     'host:s'    => \$host,
	     'dbname:s'  => \$dbname,
	     'port:n'    => \$port,
	     'module:s'  => \$module,
	     'usefile'   => \$usefile
	     );
my @est;

if( $usefile == 1 ) {
    while( <> ) {
	my ($en) = split;
	push(@est,$en);
    }
} else {
    @est = @ARGV;
}


my $locator = "$module/host=$host;port=$port;dbname=$dbname;user=$dbuser;pass=$dbpass";
print STDERR "Using $locator for todb\n";
my $db =  Bio::EnsEMBL::DBLoader->new($locator);

foreach my $est_id (@est) {
    my $trans;
    eval {
	$trans=$db->gene_Obj->get_Transcript_by_est($est_id);
    };
    if ($@) {
	print STDERR "Could not get transcript for est $est_id\n Exception: $@";
    }
    else {
	print STDERR "Got transcript ".$trans->id."\n";
    }
}
    
