#!/usr/local/bin/perl -w

use strict;
use vars qw($USER $PASS $DB $HOST $DSN);
use Bio::EnsEMBL::DBLoader;
use Bio::SeqIO;
use Bio::EnsEMBL::PerlDB::Contig;
use Bio::EnsEMBL::PerlDB::Clone;

$USER='jason';
$PASS= undef;
$DB='ensdev';
$HOST=undef;
$DSN="dbname=$DB";

$DSN .= ";user=$USER" if( defined $USER);
$DSN .= ";pass=$PASS" if( defined $PASS);
$DSN .= ";host=$HOST" if( defined $HOST);

$SIG{INT} = sub {my $sig=shift;die "exited after SIG$sig";};


my ($seqfile) = @ARGV;

if( !defined $seqfile ) { die 'cannot load because seqfile was not specified';}

my $obj = new Bio::EnsEMBL::DBLoader("Bio::EnsEMBL::DBSQL::Obj/$DSN");


my $seqio = new Bio::SeqIO(-format=>'Fasta', 
			   -file=>$seqfile);

my ($cloneid) = ($seqfile =~ /(\S+)\.\S+/);
my $count = 1;
while ( my $seq = $seqio->next_seq ) {
    my $contigid = $seq->id;
    $contigid =~ s/(\.*)\..*/$1/;
    print("Contig $contigid : contig length " . $seq->length . "\n");
       
    my $clone     = new Bio::EnsEMBL::PerlDB::Clone;
    my $contig    = new Bio::EnsEMBL::PerlDB::Contig;
    $clone->htg_phase(3);
    $clone->id($cloneid);
    $contig->id($contigid);
    $contig->internal_id($count++);
    $contig->seq($seq);    
    $clone->add_Contig($contig);
    eval { 
	my $result     = $obj->write_Clone($clone);
    };
    if( $@ ) {
	print STDERR "error was $@\n";
	last;
    }
}


