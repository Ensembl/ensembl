#!/usr/local/bin/perl

=head1 NAME

Exon duplicates test

=head1 SYNOPSIS
 
  exon_duplicates.pl

=head1 DESCRIPTION

This test checks the whole database and makes sure 
that there are no duplicated exons within a transcript

=head1 OPTIONS

    -dbtype    Database tpye (only used for TimDB)

    -host      host name for database (gets put as host= in locator)

    -port      For RDBs, what port to connect to (port= in locator)

    -dbname    For RDBs, what name to connect to (dbname= in locator)

    -dbuser    For RDBs, what username to connect as (dbuser= in locator)

    -dbpass    For RDBs, what password to use (dbpass= in locator)

    -module    Module name to load (Defaults to Bio::EnsEMBL::DBSQL::Obj)

    -help      Displays script documentation with PERLDOC

=cut

use strict;
use Bio::EnsEMBL::DBLoader;
use Bio::SeqIO;
use Getopt::Long;

my $host   = 'sol28';
my $port   = '410000';
my $dbname = 'ens100';
my $dbuser = 'ensembl';
my $dbpass = undef;
my $module = 'Bio::EnsEMBL::DBSQL::Obj';
my $help;

&GetOptions( 
	     'host:s'     => \$host,
	     'port:n'     => \$port,
	     'dbname:s'   => \$dbname,
	     'dbuser:s'   => \$dbuser,
	     'dbpass:s'   => \$dbpass,
	     'module:s'   => \$module,
	     'h|help'     => \$help
	     );

if ($help) {
    exec('perldoc', $0);
}

my $locator = "$module/host=$host;port=$port;dbname=$dbname;user=$dbuser;pass=$dbpass";
my $db =  Bio::EnsEMBL::DBLoader->new($locator);

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


