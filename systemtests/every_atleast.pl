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
	if ($clone->version == 0) {
	    $errcount++;
	    print "Error $errcount\n";
	    print "Clone $clone_id has version zero!\n";
	}
	foreach my $contig ($clone->get_all_Contigs()) {
	    print STDERR "\n        contig     ",$contig->id,"\n";
	    foreach my $gene ($contig->get_all_Genes()) {
		print STDERR "\n        gene       ",$gene->id,"\n";
		if ($gene->version == 0) {
		    $errcount++;
		    print "Error $errcount\n";
		    print "Clone:      $clone_id\n";
		    print "Contig:     ",$contig->id,"\n";
		    print "Gene ".$gene->id." has version zero!\n";
		}
		if (my $trans= $gene->each_Transcript()) {
		    foreach $trans ($gene->each_Transcript()) {
			print STDERR "\n        transcript ",$trans->id,"\n";
			my $switch = 0;
			my $start_exon_id = $trans->translation->start_exon_id;
			my $end_exon_id = $trans->translation->end_exon_id;
			if ($trans->version == 0) {
			    $errcount++;
			    print "Error $errcount\n";
			    print "Clone:      $clone_id\n";
			    print "Contig:     ",$contig->id,"\n";
			    print "Gene:       ",$gene->id,"\n";
			    print "Transcript ",$trans->id," has version 0!\n";
			}
			if (my $exon = $trans->each_Exon()) { 
			    foreach $exon ($trans->each_Exon()) {
				print STDERR "        exon       ",$exon->id,"\n";
				if ($start_exon_id eq $exon->id) {
				    $switch++;
				}
				if ($end_exon_id eq $exon->id) {
				    $switch++;
				}
				if ($exon->version == 0) {
				    $errcount++;
				    print "Error $errcount\n";
				    print "Clone:      $clone_id\n";
				    print "Contig:     ",$contig->id,"\n";
				    print "Gene:       ",$gene->id,"\n";
				    print "Transcript: ",$trans->id,"\n";
				    print "Exon ".$exon->id." has version 0!\n";
				}  
			    }
			}
			else {
			    $errcount++;
			    print "Error $errcount\n";
			    print "Clone:      $clone_id\n";
			    print "Contig:     ",$contig->id,"\n";
			    print "Gene:       ",$gene->id,"\n";
			    print "Transcript: ",$trans->id,"\n";
			    print "This transcript does not contain any exon!\n";
			}
			if ($switch != 2){
			    $errcount++;
			    print "Error $errcount\n";
			    print "Clone:      $clone_id\n";
			    print "Contig:     ",$contig->id,"\n";
			    print "Gene:       ",$gene->id,"\n";
			    print "Transcript: ",$trans->id,"\n";
			    print "This transcript does not have valid start and end exons for the translation!\n";
			}
		    }
		}
		else {
		    $errcount++;
		    print "Error $errcount\n";
		    print "Clone:      $clone_id\n";
		    print "Contig:     ",$contig->id,"\n";
		    print "Gene:       ",$gene->id,"\n";
		    print "This gene does not contain a transcript!\n";
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



