## Bioperl Test Harness Script for Modules
##
# CVS Version
# $Id$


# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

#-----------------------------------------------------------------------
## perl test harness expects the following output syntax only!
## 1..3
## ok 1  [not ok 1 (if test fails)]
## 2..3
## ok 2  [not ok 2 (if test fails)]
## 3..3
## ok 3  [not ok 3 (if test fails)]
##
## etc. etc. etc. (continue on for each tested function in the .t file)
#-----------------------------------------------------------------------


## We start with some black magic to print on failure.
BEGIN { $| = 1; print "1..2\n"; 
	use vars qw($loaded); }

END {print "not ok 1\n" unless $loaded;}

use Bio::EnsEMBL::PerlDB::Obj;

use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Exon;

use Bio::SeqIO;

use Bio::Seq;
use FileHandle;

$loaded = 1;
print "ok 1\n";    # 1st test passes.


$contig = Bio::EnsEMBL::PerlDB::Contig->new();
    
# $sf is a Bio::SeqFeatureI type object. $seq is a Bio::Seq object


$seq = Bio::Seq->new( -id => 'Contig-1' , -seq => 'ATGGCGGATGTTTATGTGGGTGGCCCGGGG' );
 
#$contig->add_SeqFeature($sf);
$contig->offset(1);
$contig->orientation(1);
$contig->order(1);
$contig->seq($seq); 
$contig->id('Contig-1');

$seq2 = Bio::Seq->new( -id => 'Contig-2' , -seq => 'TCAGAAATTTGGGTGTTTTGGCCCTGGTGGTTTGGGTTT' );
$contig2 = Bio::EnsEMBL::PerlDB::Contig->new();
$contig2->offset(900);
$contig2->orientation(1);
$contig2->order(2);
$contig2->id('Contig-2');
$contig2->seq($seq2);

$obj = Bio::EnsEMBL::PerlDB::Obj->new();

$gene = new Bio::EnsEMBL::Gene;
$tr   = new Bio::EnsEMBL::Transcript;
$ex1   = new Bio::EnsEMBL::Exon;
$ex2   = new Bio::EnsEMBL::Exon;

$ex1->start(8);
$ex1->end(13);
$ex1->phase(0);
$ex1->strand(1);
$ex1->attach_seq($seq);
$ex1->contig_id('Contig-1');
$ex1->clone_id('test-clone');
$ex1->created(time());
$ex1->modified(time());
$ex1->id('exon-id-1');

$ex2->start(18);
$ex2->end(23);
$ex2->phase(0);
$ex2->strand(1);
$ex2->attach_seq($seq2);
$ex2->contig_id('Contig-2');
$ex2->clone_id('test-clone');
$ex2->created(time());
$ex2->modified(time());
$ex2->id('exon-id-2');

$gene->id('gene-id');
$tr->id('transcript-id');

$tr->add_Exon($ex1);
$tr->add_Exon($ex2);

$trans = Bio::EnsEMBL::Translation->new();
$trans->start_exon_id('exon-id-1');
$trans->end_exon_id('exon-id-2');
$trans->start(8);
$trans->end(23);
$tr->translation($trans);
$tr->id('peptide-id');
$gene->add_Transcript($tr);

$contig->add_Gene($gene);

$obj->write_Contig($contig);
$obj->write_Contig($contig2);


$clone = Bio::EnsEMBL::PerlDB::Clone->new();
$clone->id("test-clone");

$clone->add_Contig($contig);
$clone->add_Contig($contig2);


print "ok 2\n";
exit(0);
$fh = new FileHandle;
$fh->open('>out.embl');

$as = $clone->get_AnnSeq();
$asio = Bio::SeqIO->new( '-format' => 'EMBL' , -fh => $fh) ;
$asio->_post_sort(\&sort_FTHelper_EnsEMBL);

$asio->write_seq($as);

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











