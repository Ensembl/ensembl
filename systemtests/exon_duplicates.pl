#!/usr/local/bin/perl

=head1 NAME

Exon duplicates

=head1 SYNOPSIS
 
  exon_duplicates.pl

=head1 DESCRIPTION

This test checks the whole database and makes sure 
that there are no duplicated exons within a transcript

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
		foreach my $trans ($gene->each_Transcript()) {
		    print STDERR "\n        transcript ",$trans->id,"\n";
		    my $exon, my $i, my $j;
		    my $switch = 0;
		    my @starts=();
		    my @ends=();
		    my @ids=();
		    foreach $exon ($trans->each_Exon()) {
			print STDERR "        exon       ",my $id=$exon->id;
			print STDERR " (start: ",my $start = $exon->start()," end: ",my $end = $exon->end(),")\n";
			push (@starts, $start);
			push (@ends, $end);
			push (@ids, $id);
		    }
		    for ($i=0; $i<=@starts;$i++) {
			for ($j=$i+1; $j<=@starts;$j++) {
			    if ($starts[$i] == $starts[$j] && $ends[$i]==$ends[$j]) {
				$errcount++;
				print "Error $errcount\n";
				print "Clone:      $clone_id\n";
				print "Contig:     ",$contig->id,"\n";
				print "Gene:       ",$gene->id,"\n";
				print "Transcript: ",$trans->id,"\n";
				print "Exon $ids[$i] has the same start and and as exon $ids[$j] \n";
			    }
			}
		    } 
		}
	    }
	}
    };

    if( $@ ) {
	print "Unable to process $clone_id due to \n$@\n";
    }
}

if ($errcount>0) {
    print STDERR "\nFound $errcount  duplicated exons\n";
}


