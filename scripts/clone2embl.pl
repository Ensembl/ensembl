#!/usr/local/bin/perl

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
# this doesn't have genes (finished)
#my $clone  = 'dJ1156N12';
# this does have genes (finished)
my $clone  = 'dJ271M21';
# this does have genes (unfinished)
# my $clone = '217N14';

&GetOptions( 'dbtype:s' => \$dbtype,
	     'host:s'   => \$host,
	     'port:n'   => \$port,
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

$as->seq->desc("Reannotated Clone via EnsEMBL");
my $comment = Bio::Annotation::Comment->new();

$comment->text("This clone was reannotated via the EnsEMBL system. Please visit the EnsEMBL web site, http://ensembl.ebi.ac.uk for more information");

$as->annotation->add_Comment($comment);

my $emblout = Bio::AnnSeqIO->new( -format => 'EMBL', -fh => \*STDOUT);

$emblout->write_annseq($as);


