#!/usr/local/bin/perl

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
	print STDERR "\nDumping clone      $clone_id\n";
    }
    
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
			    print ERROR "Error $errcount\n";
			    print ERROR "Clone:      $clone_id\n";
			    print ERROR "Contig:     ",$contig->id,"\n";
			    print ERROR "Gene:       ",$gene->id,"\n";
			    print ERROR "Transcript: ",$trans->id,"\n";
			    print ERROR "This transcript does not contain any exon!\n";
			}   
		    }
		}
		else {
		    $errcount++;
		    print ERROR "Error $errcount\n";
		    print ERROR "Clone:      $clone_id\n";
		    print ERROR "Contig:     ",$contig->id,"\n";
		    print ERROR "Gene:       ",$gene->id,"\n";
		    print ERROR "This gene does not contain a transcript!\n";
		}
	    }
	}
    };
    if( $@ ) {
	print ERROR "Unable to process $clone_id due to \n$@\n";
    }
    
    if ($errcount>0) {
	print STDERR "\nFound $errcount empty genes/transcripts (see error.log for details)\n";
    }
}
close ERROR;

