#!/usr/local/bin/perl

=head1 NAME

DNA test

=head1 SYNOPSIS
 
  dna_test.pl

=head1 DESCRIPTION

This test checks the whole database and makes sure 
that all clones contain dna.

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
    print STDERR "\nDumping clone $clone_id:\n";

    eval {
	my $clone = $db->get_Clone($clone_id);
	my @contig_id = $clone->get_all_Contigs();
	foreach my $contig (@contig_id) {
	    print STDERR "       contig ".$contig->id."\n";
	    my $contig_seq = $contig->seq->seq();
	    if ( $contig_seq =~ /[^A,T,G,C,N,R,Y]/ || $contig_seq eq "") {
		$errcount++;
		$contig_seq =~ s/[A,T,G,C,N]//g;
		print "Error $errcount\n";
		print "Clone:   $clone_id\n";
		print "Contig:  ",$contig->id,"\n";
		print "Error:\n";
		if ($contig_seq eq "") {
		    print "no sequence present in this contig!\n";
		}
		else {
		    print "non-DNA sequence found:\"$contig_seq\"\n\n";
		}
		next;
	    }
	}

    };

    if( $@ ) {
	print "Unable to process $clone_id due to \n$@\n";
    }
}

if ($errcount>0) {
	print STDERR "\nFound $errcount contig(s) not containing DNA\n";
    }

