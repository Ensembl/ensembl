#!/usr/local/bin/perl

=head1 NAME

Transcript strand

=head1 SYNOPSIS
 
  transcript_strand.pl

=head1 DESCRIPTION

This test checks the whole database and makes sure 
that exons of transcripts on the same internal contig 
are on the same strand. 

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

my $host   = 'localhost';
my $port   = '410000';
my $dbname = 'ensdev';
my $dbuser = 'ensembl';
my $dbpass = undef;
my $module = 'Bio::EnsEMBL::DBOLD::Obj';
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
				print "Error $errcount\n";
				print "Clone:      $clone_id\n";
				print "Contig:     ",$contig->id,"\n";
				print "Gene:       ",$gene->id,"\n";
				print "Transcript: ",$trans->id,"\n";
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
	print "Unable to process $clone_id due to \n$@\n";
    }
}

if ($errcount>0) {
	print STDERR "\nFound $errcount transcript(s) with exons in wrong strands(see error.log for details)\n";
    }

