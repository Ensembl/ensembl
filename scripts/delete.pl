#!/usr/local/bin/perl

use strict;

use Getopt::Long;
use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DBLoader;

# global defaults
my $host = 'localhost';
my $dbuser = 'ensdev';
my $dbname = 'ensembl';
my $port = '410000';
my $dbpass = undef;
my $do_gene = 1;
my $usefile = 0;
my $module = 'Bio::EnsEMBL::DBSQL::Obj';

&GetOptions( 
	     'dbuser:s'  => \$dbuser,
	     'usefile:s' => \$usefile,
	     'dbpass:s'  => \$dbpass,
	     'host:s'    => \$host,
	     'dbname:s'  => \$dbname,
	     );

my $locator = "$module/host=$host;port=$port;dbname=$dbname;user=$dbuser;pass=$dbpass";
my $db =  Bio::EnsEMBL::DBLoader->new($locator);


my @clones;
if ($usefile) {
    print STDERR "Using ".$usefile." as a list of clones to delete\n";
    open(IN,"<$usefile");
    while(<IN>){
	if(/^(\S+)/){
	    push(@clones,$1);
	}
    }
}
else {
    @clones=@ARGV;
}

foreach my $clone_id ( @clones) {
    print STDERR "Deleting $clone_id\n";
    eval {
	my $clone = $db->get_Clone($clone_id);
	
	if( $do_gene == 1 ) {
	    my @genes = $clone->get_all_Genes();
	    
	    foreach my $gene ( @genes ) {
		$db->delete_Gene($gene->id());
	    }
	}
	
	$db->delete_Clone($clone_id);
    };
    if ($@) {
	print STDERR "Could not delete clone $clone_id:\n$@\n";
    }
}


	

