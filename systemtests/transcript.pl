#!usr/local/bin/perl
=head1 NAME

Every_atleast

=head1 SYNOPSIS
 
  every_atleast.pl

=head1 DESCRIPTION

This test checks the whole database and makes sure 
that every gene has at least one transcript, and 
every transcript at least one exon. Also, it checks 
that every transcript with a translation has a 
correct start and end exons.

=head1 OPTIONS

    -dbtype    Database tpye (only used for TimDB)

    -host      host name for database (gets put as host= in locator)

    -port      For RDBs, what port to connect to (port= in locator)

    -dbname    For RDBs, what name to connect to (dbname= in locator)

    -dbuser    For RDBs, what username to connect as (dbuser= in locator)

    -dbpass    For RDBs, what password to use (dbpass= in locator)

    -module    Module name to load (Defaults to Bio::EnsEMBL::DBSQL::Obj)

    -help      Displays script documentation with PERLDOC

    -usefile   read in on stdin a list of clones, one clone per line

    -getall    get all clones
=cut

use strict;
use Bio::EnsEMBL::DBLoader;
use Bio::SeqIO;
use Getopt::Long;

my $dbtype = 'rdb';
my $host   = 'sol28';
my $port   = '410000';
my $dbname = 'ens100';
my $dbuser = 'ensembl';
my $dbpass = undef;
my $module = 'Bio::EnsEMBL::DBSQL::Obj';
my $help;
my $usefile = 0;
my $getall = 0;

&GetOptions( 
	     'dbtype:s'   => \$dbtype,
	     'host:s'     => \$host,
	     'port:n'     => \$port,
	     'dbname:s'   => \$dbname,
	     'dbuser:s'   => \$dbuser,
	     'dbpass:s'   => \$dbpass,
	     'module:s'   => \$module,
	     'usefile'    => \$usefile,
	     'getall'     => \$getall,
	     'h|help'     => \$help
	     );


if ($help) {
    exec('perldoc', $0);
}

my $locator = "$module/host=$host;port=$port;dbname=$dbname;user=$dbuser;pass=$dbpass";
my $db =  Bio::EnsEMBL::DBLoader->new($locator);
my $seqio;
my $errcount = 0;
my @clone_id;

if( $usefile == 1 ) {
    while( <> ) {
	my ($en) = split;
	push(@clone_id,$en);
    }
} 

elsif ($getall == 1) {
    @clone_id = $db->get_all_Clone_id();
}

else {
    @clone_id = @ARGV;
} 

foreach my $clone_id ( @clone_id ) {
    print STDERR "\nDumping clone      $clone_id\n";
    eval {
	my $clone = $db->get_Clone($clone_id);
	foreach my $contig ($clone->get_all_Contigs()) {
	    print STDERR "\n        contig     ",$contig->id,"\n";
	    foreach my $gene ($contig->get_all_Genes()) {
		print STDERR "\n        gene       ",$gene->id,"\n";
		foreach my $trans ($gene->each_Transcript()) {
		    print STDERR "\n        transcript ",$trans->id,"\n";
		    my $trans_seq = $trans->dna_seq();
		    my $err;
		    if ($trans_seq eq "") {
			$err = "no sequence present in this contig!\n";
		    }
		    $trans_seq =~ s/[A,T,G,C,N,Y,R]//g;
		    print "$trans_seq\n";
		    if ($trans_seq ne "") {
			$errcount++;
			print "Error $errcount\n";
			print "Clone:     $clone_id\n";
			print "Contig:    ",$contig->id,"\n";
			print "Gene:      ",$gene->id,"\n";
			print "Transcript ",$trans->id,"\n";
			print "Error:\n";
			if ($err) {
			    print $err;
			}
			else {
			    print "non-DNA sequence found:\"$trans_seq\"\n\n";
			}
		    }
		}
	    }
	}
    };
    if( $@ ) {
	print "Unable to process $clone_id due to \n$@\n";
    }
    
    if ($errcount>0) {
	print STDERR "\nFound $errcount errors\n";
    }
}



