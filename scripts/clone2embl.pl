#!/usr/local/bin/perl

use strict;

use Bio::EnsEMBL::AceDB::Obj;
use Bio::EnsEMBL::DB::Obj;
use Bio::EnsEMBL::TimDB::Obj;
use Bio::AnnSeqIO;

use Getopt::Long;

my $dbtype = 'rdb';
my $host   = 'croc';
my $host1   = 'humsrv1';
my $port   = '410000';
my $clone  = 'dJ1156N12';

&GetOptions( 'dbtype:s' => \$dbtype,
	     'host:s'   => \$host,
	     'port:n'   => \$port,
	     );

my $db;

if( $dbtype =~ 'ace' ) {
    $db = Bio::EnsEMBL::AceDB::Obj->new( -host => $host1, -port => $port);
} elsif ( $dbtype =~ 'rdb' ) {
    $db = Bio::EnsEMBL::DB::Obj->new( -user => 'root', -db => 'pog' , -host => $host );
} elsif ( $dbtype =~ 'timdb' ) {
    $db = Bio::EnsEMBL::TimDB::Obj->new();
} else {
    die("$dbtype is not a good type (should be ace, rdb or timdb)");
}


my $clone_id = shift;
$clone_id=$clone unless $clone_id;


my $clone = $db->get_Clone($clone_id);
my $as = $clone->get_AnnSeq();


$as->seq->desc("Reannotated Clone via EnsEMBL");
my $comment = Bio::Annotation::Comment->new();

$comment->text("This clone was reannotated via the EnsEMBL system. Please visit the EnsEMBL web site, http://ensembl.ebi.ac.uk for more information");

$as->annotation->add_Comment($comment);


my $emblout = Bio::AnnSeqIO->new( -format => 'EMBL', -fh => \*STDOUT);

$emblout->write_annseq($as);


