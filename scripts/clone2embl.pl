#!/usr/local/bin/perl

=head1 NAME - clone2embl

 provides flat file formats from EnsEMBL databases

=head1 SYNOPSIS - 

    clone2embl dJ271M21

    clone2embl -gff dJ271M21
   
    clone2embl -dbtype ace dJ271M21

    clone2embl -dbtype rdb -host mysql.server.somewhere dJ271M21

=cut


use strict;

use Bio::EnsEMBL::AceDB::Obj;
use Bio::EnsEMBL::DB::Obj;
use Bio::EnsEMBL::TimDB::Obj;
use Bio::AnnSeqIO;

use Getopt::Long;

my $dbtype = 'rdb';
my $host;
my $host1  = 'croc';
my $host2  = 'humsrv1';
my $port   = '410000';
my $gff;
my $nodna = 0;

# this doesn't have genes (finished)
#my $clone  = 'dJ1156N12';
# this does have genes (finished)
my $clone  = 'dJ271M21';
# this does have genes (unfinished)
# my $clone = '217N14';

&GetOptions( 'dbtype:s' => \$dbtype,
	     'host:s'   => \$host,
	     'port:n'   => \$port,
	     'gff'      => \$gff,
	     'nodna'      => \$nodna,
	     );

my $db;

my $clone_id = shift;
$clone_id=$clone unless $clone_id;

if( $dbtype =~ 'ace' ) {
    $host=$host2 unless $host;
    $db = Bio::EnsEMBL::AceDB::Obj->new( -host => $host, -port => $port);
} elsif ( $dbtype =~ 'rdb' ) {
    $host=$host1 unless $host;
    $db = Bio::EnsEMBL::DB::Obj->new( -user => 'root', -db => 'pog' , -host => $host );
} elsif ( $dbtype =~ 'timdb' ) {
    $db = Bio::EnsEMBL::TimDB::Obj->new($clone_id);
} else {
    die("$dbtype is not a good type (should be ace, rdb or timdb)");
}

my $clone = $db->get_Clone($clone_id);
my $as = $clone->get_AnnSeq();

# choose output mode

if($gff){
    
    # only works with one contig for now
    if(scalar($clone->get_all_Contigs)!=1){
	$clone->throw("More than one contig in clone ".$clone->id()."\n");
    }
    my ($contig)=$clone->get_all_Contigs;
    my $gff=$contig->gff;
    $gff->dump;

}else{
    $as->seq->desc("Reannotated Clone via EnsEMBL");
    my $comment = Bio::Annotation::Comment->new();

    $comment->text("This clone was reannotated via the EnsEMBL system. Please visit the EnsEMBL web site, http://ensembl.ebi.ac.uk for more information");

    $as->annotation->add_Comment($comment);

    my $emblout = Bio::AnnSeqIO->new( -format => 'EMBL', -fh => \*STDOUT);
    $emblout->_post_sort(\&sort_FTHelper_EnsEMBL);
    if( $nodna == 1 ) {
	$emblout->_show_dna(0);
    }

    

    $emblout->write_annseq($as);
}

sub sort_FTHelper_EnsEMBL {
    my $a = shift;
    my $b = shift;

    if( $a->key eq $b->key ) {
	return ($a->loc cmp $b->loc);
    }

    if( $a->key eq 'CDS' ) {
	return -1;
    }

    return 1;
}

