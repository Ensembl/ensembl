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

=cut

use strict;
use Bio::EnsEMBL::DBSQL::Obj;
use Bio::SeqIO;
use Getopt::Long;

my $tdbtype = 'rdb';
my $thost   = 'sol28';
my $tport   = '410000';
my $tdbname = 'ensdev';
my $user    = 'ensembl';

&GetOptions( 
	     'dbtype:s'   => \$tdbtype,
	     'host:s'     => \$thost,
	     'port:n'     => \$tport,
	     'user:s'     => \$user,
	     'dbname:s'   => \$tdbname,
	     );

my $db = Bio::EnsEMBL::DBSQL::Obj->new( -user => $user, -db => $tdbname , -host => $thost );
my @clone_id = $db->get_all_Clone_id();
my $seqio;
my $errcount = 0;

foreach my $clone_id ( @clone_id ) {
    print STDERR "\nDumping clone      $clone_id\n";
    eval {
	my $clone = $db->get_Clone($clone_id);
	foreach my $contig ($clone->get_all_Contigs()) {
	    print STDERR "\n        contig     ",$contig->id,"\n";
	    foreach my $gene ($contig->get_all_Genes()) {
		print STDERR "\n        gene       ",$gene->id,"\n";
		if (my $trans= $gene->each_Transcript()) {
		    foreach $trans ($gene->each_Transcript()) {
			print STDERR "\n        transcript ",$trans->id,"\n";
			if (my $exon = $trans->each_Exon()) { 
			    foreach $exon ($trans->each_Exon()) {
				print STDERR "        exon       ",$exon->id,"\n";
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
	print STDERR "\nFound $errcount empty genes/transcripts\n";
    }
}


