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
my $format  = 'pep';
my $usefile = 0;
my $verbose = 1;
my $noacc   = 0;
my $test    = 0;
my $user    = 'ensembl';
my $logerror = 'error.log';

&GetOptions( 
	     'host:s'     => \$thost,
	     'port:n'     => \$tport,
	     'user:s'     => \$user,
	     'usefile'    => \$usefile,
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

my @gene_id = $db->get_all_Gene_id();

my $seqio;

if( $format eq 'pep' ) {
    $seqio = Bio::SeqIO->new('-format' => 'Fasta' , -fh => \*STDOUT ) ;
}

#@gene_id = ('dummy_gene');
my $errcount = 0;

foreach my $gene_id ( @gene_id ) {
    if( $verbose ) {
	print STDERR "Dumping $gene_id\n";
    }

    eval {
	my $gene = $db->get_Gene($gene_id);
	foreach my $trans ( $gene->each_Transcript ) {

	    # get out first exon. Tag it to clone and gene on this basis
	    my @exon = $trans->each_Exon;
	    my $fe = $exon[0];
	    my $tseq = $trans->translate()->seq();
	    if ( $tseq =~ /\*/ ) {
		$errcount++;
		print ERROR "Stop codons found in the translation of:\n";
		print ERROR "Clone:          ". $fe->clone_id ."\n";
		print ERROR "Gene:           $gene_id\n";
		print ERROR "Transcript:     ",$trans->id,"\n";
		next;
	    }
	}

    };
    if ($errcount>0) {
	print STDERR "\nFound $errcount transcript(s) translating into a stop codon (see error.log for details)\n";
    }
    if( $@ ) {
	print ERROR "Unable to process $gene_id due to \n$@\n";
    }
}

close ERROR;
