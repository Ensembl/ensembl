#!/usr/local/bin/perl

=head1 NAME

Data integrity test 1

=head1 SYNOPSIS
 
  test1.pl

=head1 DESCRIPTION

This test checks the whole database and makes sure 
that the translations of all genes do not contain 
any stop codons.

=head1 OPTIONS

   -getall 


=cut
use strict;
use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::TimDB::Obj;
use Bio::SeqIO;

use Getopt::Long;

my $thost   = 'sol28';
my $tport   = '410000';
my $tdbname = 'ensdev';
my $user    = 'ensembl';

&GetOptions( 
	     'host:s'     => \$thost,
	     'port:n'     => \$tport,
	     'user:s'     => \$user,
	     'dbname:s'   => \$tdbname,
	     );
my $db = Bio::EnsEMBL::DBSQL::Obj->new( -user => $user, -db => $tdbname , -host => $thost );
my @gene_id = $db->get_all_Gene_id();
my $seqio;
my $errcount = 0;

foreach my $gene_id ( @gene_id ) {
    print STDERR "Dumping $gene_id\n";

    eval {
	my $gene = $db->get_Gene($gene_id);
	foreach my $trans ( $gene->each_Transcript ) {
	    my @exon = $trans->each_Exon;
	    my $fe = $exon[0];
	    my $tseq = $trans->translate()->seq();
	    if ( $tseq =~ /\*/ ) {
		$errcount++;
		print "Stop codons found in the translation of:\n";
		print "Clone:          ". $fe->clone_id ."\n";
		print "Gene:           $gene_id\n";
		print "Transcript:     ",$trans->id,"\n";
		next;
	    }
	}

    };
    
    if ($errcount>0) {
	print STDERR "\nFound $errcount transcript(s) translating into a stop codon\n";
    }

    if( $@ ) {
	print "Unable to process $gene_id due to \n$@\n";
    }
}

