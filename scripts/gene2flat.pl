#!/usr/local/bin/perl

=head1 NAME

gene2flat

=head1 SYNOPSIS
 
  gene2flat ENSG00000012

=head1 DESCRIPTION

gene2flat produces a number of flat file outputs of the genes,
in particular the protein translation

=head1 OPTIONS

   -getall 


=cut
use strict;
use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::TimDB::Obj;
use Bio::SeqIO;

use Getopt::Long;

my $tdbtype = 'rdb';
my $thost   = 'sol28';
my $tport   = '410000';
my $tdbname = 'ensdev';
my $format  = 'transcript';
my $usefile = 0;
my $getall  = 1;
my $verbose = 0;
my $noacc   = 0;
my $test    = 0;
my $user    = 'ensembl';
my $logerror = undef;

&GetOptions( 
	     'dbtype:s'   => \$tdbtype,
	     'host:s'     => \$thost,
	     'port:n'     => \$tport,
	     'user:s'     => \$user,
	     'usefile'    => \$usefile,
	     'dbname:s'   => \$tdbname,
	     'format:s'   => \$format,
	     'getall'     => \$getall,
	     'verbose'    => \$verbose,
	     'test'       => \$test,
	     'noacc'      => \$noacc,
	     'logerror:s'   => \$logerror,
	     );
my $db;

if( defined $logerror ) {
    open(ERROR,">$logerror") || die "Could not open $logerror $!";
} else {
    open(ERROR,">&STDERR") || die "Could not dup STDERR";
}


if( $tdbtype =~ 'ace' ) {
    $db = Bio::EnsEMBL::AceDB::Obj->new( -host => $thost, -port => $tport);
} elsif ( $tdbtype =~ 'rdb' ) {
    $db = Bio::EnsEMBL::DBSQL::Obj->new( -user => $user, -db => $tdbname , -host => $thost );
} elsif ( $tdbtype =~ 'timdb' ) {
    $db = Bio::EnsEMBL::TimDB::Obj->new('',$noacc,$test);
} else {
    die("$tdbtype is not a good type (should be ace, rdb, timdb)");
}

my @gene_id;

if( $usefile ) {
    while( <> ) {
	my ($g) = split;
	push(@gene_id,$g);
    }
} elsif ( $getall == 1 ) {
    @gene_id = $db->get_all_Gene_id();
} else {
    @gene_id = @ARGV;
}

my $seqio;

if( $format eq 'pep' || $format eq 'transcript' ) {
    $seqio = Bio::SeqIO->new('-format' => 'Fasta' , -fh => \*STDOUT ) ;
}

foreach my $gene_id ( @gene_id ) {
    if( $verbose ) {
	print STDERR "Dumping $gene_id\n";
    }

    eval {

	my $gene = $db->get_Gene($gene_id);

	if( $format eq 'pep' ) {
	    foreach my $trans ( $gene->each_Transcript ) {
		# get out first exon. Tag it to clone and gene on this basis
		my @exon = $trans->each_Exon;
		my $fe = $exon[0];
		my $tseq = $trans->translate();
		if ( $tseq->seq =~ /\*/ ) {
		    print ERROR "translation has stop codons. Skipping! (in clone". $fe->clone_id .")\n";
		    next;
		}
		$tseq->desc("Gene:$gene_id Clone:".$fe->clone_id);
		$seqio->write_seq($tseq);
	    }
	} elsif ( $format eq 'dump' ) {
	    foreach my $trans ( $gene->each_Transcript ) {
		print "Transcript ",$trans->id,"\n";
		foreach my $exon ( $trans->each_Exon ) {
		    print "  Exon ",$exon->id," ",$exon->contig_id,":",$exon->start,"-",$exon->end,".",$exon->strand,"\n";
		    my $seq = $exon->seq();
		    my $str = $seq->str();
		    print "    Start phase ",$exon->phase,"[",substr($str,0,10),"] End phase ",$exon->end_phase," [",substr($str,-10),"]\n";
		}
	    }

	} 
	elsif ($format eq 'transcript') {
	    foreach my $trans ( $gene->each_Transcript ) {
		my $seq = $trans->dna_seq();
		$seq->id($trans->id);
		my @exon = $trans->each_Exon;
		my $fe = $exon[0];
		$seq->desc("Gene:$gene_id Clone:".$fe->clone_id);
		$seqio->write_seq($seq);
	    }
	}
	else {
	    die "No valid format!";
	}
    };

    if( $@ ) {
	print ERROR "Unable to process $gene_id due to \n$@\n";
    }
}





