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
use Bio::EnsEMBL::DBLoader;
use Bio::SeqIO;
use Getopt::Long;
use vars qw(@ISA);

@ISA = qw(Bio::Root::Object);

my $dbtype = 'rdb';
my $host   = 'obi-wan';
my $port   = '410000';
my $dbname = 'ensembl';
my $dbuser = 'ensro';
my $dbpass = undef;
my $module = 'Bio::EnsEMBL::DBSQL::Obj';
my $help;

&GetOptions( 
	     'dbtype:s'   => \$dbtype,
	     'host:s'     => \$host,
	     'port:n'     => \$port,
	     'dbname:s'   => \$dbname,
	     'dbuser:s'   => \$dbuser,
	     'dbpass:s'   => \$dbpass,
	     'module:s'   => \$module,
	     'h|help'     => \$help,
	     );


if ($help) {
    exec('perldoc', $0);
}

print STDERR "\nConnecting to $dbname database...\n";
my $recipient_locator = "$module/host=$host;port=$port;dbname=$dbname;user=$dbuser;pass=$dbpass";
my $db =  Bio::EnsEMBL::DBLoader->new($recipient_locator);

my @clone_id = $db->get_all_Clone_id();
my $seqio;
my $errcount = 0;
my $err;

foreach my $clone_id ( @clone_id ) {
    print STDERR "\nDumping clone $clone_id:\n";

    eval {
	my $clone = $db->get_Clone($clone_id);
	my @contig_id = $clone->get_all_Contigs();
	foreach my $contig (@contig_id) {
	    print STDERR "       contig ".$contig->id."\n";
	    my $contig_seq = $contig->seq->seq();
	    
	    if ($contig_seq eq "") {
		$err = "no sequence present in this contig!\n";
	    }
	    $contig_seq =~ s/[A,T,G,C,N,Y,R]//g;
	    print "$contig_seq\n";
	    if ($contig_seq ne "") {
		$errcount++;
		print "Error $errcount\n";
		print "Clone:   $clone_id\n";
		print "Contig:  ",$contig->id,"\n";
		print "Error:\n";
		if ($err) {
		    print $err;
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

