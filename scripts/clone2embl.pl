#!/usr/local/bin/perl

use strict;

use Bio::EnsEMBL::AceDB::Obj;
use Bio::AnnSeqIO;

my $db = Bio::EnsEMBL::AceDB::Obj->new( -host => 'humsrv1', -port => '410000');

$db->_exon_id_start("HE00000040");



my $clone_id = shift;

my $clone = $db->get_Clone($clone_id);

my $as = $clone->get_AnnSeq();


$as->seq->desc("Reannotated Clone via EnsEMBL");
my $comment = Bio::Annotation::Comment->new();

$comment->text("This clone was reannotated via the EnsEMBL system. Please visit the EnsEMBL web site, http://ensembl.ebi.ac.uk for more information");

$as->annotation->add_Comment($comment);


my $emblout = Bio::AnnSeqIO->new( -format => 'EMBL', -fh => \*STDOUT);

$emblout->write_annseq($as);


