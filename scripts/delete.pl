#!/usr/local/bin/perl

use strict;

use Getopt::Long;
use Bio::EnsEMBL::DBSQL::Obj;

# global defaults
my $host = 'localhost';
my $dbuser = 'ensadmin';
my $dbname = 'ensembl';
my $dbpass = undef;
my $do_gene = 1;


&GetOptions( 
	     'dbuser:s'  => \$dbuser,
	     'dbpass:s'  => \$dbpass,
	     'host:s'    => \$host,
	     'dbname:s'  => \$dbname,
	     );

my $db = Bio::EnsEMBL::DBSQL::Obj->new( -user => $dbuser, -dbname => $dbname , 
					-host => $host, -password => $dbpass );



foreach my $clone_id ( @ARGV) {
    print STDERR "Deleting $clone_id\n";
    my $clone = $db->get_Clone($clone_id);

    if( $do_gene == 1 ) {
	my @genes = $clone->get_all_Genes();

	foreach my $gene ( @genes ) {
	    $db->delete_Gene($gene->id());
	}
    }

    $db->delete_Clone($clone_id);
    
}

	

