#!/usr/local/bin/perl

use strict;

use Getopt::Long;
use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DBLoader;

# global defaults
my $host = 'localhost';
my $dbuser = 'ensadmin';
my $dbname = 'ensembl';
my $dbpass = undef;
my $module = 'Bio::EnsEMBL::DBSQL::Obj';
my $port = '410000';
my $do_gene = 1;
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
my @clone;

if( $usefile == 1 ) {
    while( <> ) {
	my ($en) = split;
	push(@clone,$en);
    }
} else {
    @clone = @ARGV;
}


my $locator = "$module/host=$host;port=$port;dbname=$dbname;user=$dbuser;pass=$dbpass";
print STDERR "Using $locator for todb\n";
my $db =  Bio::EnsEMBL::DBLoader->new($locator);

foreach my $clone_id ( @clone) {
    print STDERR "Deleting $clone_id\n";
    my $clone;
    eval {
	$clone = $db->get_Clone($clone_id);
    };
    if ($@) {
	print STDERR "Skipping clone $clone_id, not present in db!$@\n";
	next;
    }
    if( $do_gene == 1 ) {
	my @genes = $clone->get_all_Genes();

	foreach my $gene ( @genes ) {
	    print STDERR "Deleting gene ".$gene->id."\n";
	    $db->delete_Gene($gene->id());
	}
    }
    $db->delete_Clone($clone_id);
    
}

	

