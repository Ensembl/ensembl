#!/usr/local/bin/perl

=head1 NAME

Transcript strand

=head1 SYNOPSIS
 
  transcript_strand.pl

=head1 DESCRIPTION

This test checks the whole database and makes sure 
that exons of transcripts on the same internal contig 
are on the same strand. 

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
	 TRANS: foreach my $trans ($gene->each_Transcript()) {
		    print STDERR "\n        transcript ",$trans->id,"\n";
		    my $exon;
		    my $switch = 0;
		    my $prev_strand;
		    my $strand;
		    my $prev_contid;
		    my $contid;
		    foreach $exon ($trans->each_Exon()) {
			print STDERR "        exon       ",$exon->id;
			print STDERR "        (strand ",$exon->strand,")\n";
			if ($switch == 0) {
			    $prev_strand = $exon->strand;
			    $prev_contid = $exon->contig_id;
			}
			else {
			    $strand = $exon->strand;
			    $contid = $exon->contig_id;
			    if ($strand != $prev_strand && $contid eq $prev_contid) {
				$errcount++;
				print ERROR "Error $errcount\n";
				print ERROR "Clone:      $clone_id\n";
				print ERROR "Contig:     ",$contig->id,"\n";
				print ERROR "Gene:       ",$gene->id,"\n";
				print ERROR "Transcript: ",$trans->id,"\n";
				last TRANS;
			    }
			    $prev_strand = $strand;
			    $prev_contid = $contid;
			}
			$switch =1;
		    }
		}
	    }
	}
    };
    if( $@ ) {
	print ERROR "Unable to process $clone_id due to \n$@\n";
    }
}
if ($errcount>0) {
	print STDERR "\nFound $errcount transcript(s) with exons in wrong strands(see error.log for details)\n";
    }

close ERROR;

