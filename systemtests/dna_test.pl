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
my $format  = 'pep';
my $verbose = 1;
my $noacc   = 0;
my $test    = 0;
my $user    = 'ensembl';
my $logerror = 'error.log';

&GetOptions( 
	     'dbtype:s'   => \$tdbtype,
	     'host:s'     => \$thost,
	     'port:n'     => \$tport,
	     'user:s'     => \$user,
	     'dbname:s'   => \$tdbname,
	     'format:s'   => \$format,
	     'verbose'    => \$verbose,
	     'test'       => \$test,
	     'noacc'      => \$noacc,
	     'logerror:s'   => \$logerror,
	     );
my $db;


open(ERROR,">$logerror") || die "Could not open $logerror $!";

$db = Bio::EnsEMBL::DBSQL::Obj->new( -user => $user, -db => $tdbname , -host => $thost );

my @clone_id = $db->get_all_Clone_id();

my $seqio;

if( $format eq 'pep' ) {
    $seqio = Bio::SeqIO->new('-format' => 'Fasta' , -fh => \*STDOUT ) ;
}

#@clone_id = ('dummy_clone');

my $errcount = 0;

foreach my $clone_id ( @clone_id ) {
    if( $verbose ) {
	print STDERR "\nDumping clone $clone_id:\n";
    }

    eval {
	my $clone = $db->get_Clone($clone_id);
	my @contig_id = $clone->get_all_Contigs();
	foreach my $contig (@contig_id) {
	    # get out first exon. Tag it to clone and clone on this basis
	    print STDERR "       contig ".$contig->id."\n";
	    my $contig_seq = $contig->seq->seq();
	    if ( $contig_seq =~ /[^A,T,G,C,N,R,Y]/ || $contig_seq eq "") {
		$errcount++;
		$contig_seq =~ s/[A,T,G,C,N]//g;
		print ERROR "Error $errcount\n";
		print ERROR "Clone:   $clone_id\n";
		print ERROR "Contig:  ",$contig->id,"\n";
		print ERROR "Error:\n";
		if ($contig_seq eq "") {
		    print ERROR "no sequence present in this contig!\n";
		}
		else {
		    print ERROR "non-DNA sequence found:\"$contig_seq\"\n\n";
		}
		next;
	    }
	}

    };
    if( $@ ) {
	print ERROR "Unable to process $clone_id due to \n$@\n";
    }
}
if ($errcount>0) {
	print STDERR "\nFound $errcount transcript(s) not containing DNA (see error.log for details)\n";
    }

close ERROR;

