#!/usr/local/bin/perl

use strict;
use Bio::EnsEMBL::DBSQL::Obj;
use Bio::SeqIO;

use Getopt::Long;

my $tdbtype = 'rdb';
my $thost   = 'croc';
my $tport   = '410000';
my $tdbname = 'ensdev';
my $format  = 'pep';
my $usefile = 0;

&GetOptions( 
	     'dbtype:s' => \$tdbtype,
	     'host:s'   => \$thost,
	     'port:n'   => \$tport,
	     'usefile'  => \$usefile,
	     'dbname:s' => \$tdbname,
	     'format:s'   => \$format,
	     );
my $db;

if( $tdbtype =~ 'ace' ) {
    $db = Bio::EnsEMBL::AceDB::Obj->new( -host => $thost, -port => $tport);
} elsif ( $tdbtype =~ 'rdb' ) {
    $db = Bio::EnsEMBL::DBSQL::Obj->new( -user => 'root', -db => $tdbname , -host => $thost );
} else {
    die("$tdbtype is not a good type (should be ace, rdb)");
}

my @gene_id;

if( $usefile ) {
    while( <> ) {
	my ($g) = split;
	push(@gene_id,$g);
    }
} else {
    @gene_id = @ARGV;
}

my $seqio;

if( $format eq 'pep' ) {
    $seqio = Bio::SeqIO->new('-format' => 'Fasta' , -fh => \*STDOUT ) ;
}

foreach my $gene_id ( @gene_id ) {

    eval {

	my $gene = $db->get_Gene($gene_id);

	if( $format eq 'pep' ) {
	    foreach my $trans ( $gene->each_Transcript ) {
		my $tseq = $trans->translate();
		$seqio->write_seq($tseq);
	    }
	} else {
	    die "No valid format!";
	}
    };

    if( $@ ) {
	print STDERR "Unable to process $gene_id due to \n$@\n";
    }
}
