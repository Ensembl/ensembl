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
BEGIN { $| = 1; print "1..9\n"; 
	use vars qw($loaded); }

END {print "not ok 1\n" unless $loaded;}

use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Translation;

$loaded = 1;
print "ok 1\n";    # 1st test passes.

$gene = new Bio::EnsEMBL::Gene;
print "ok 2\n";   
$tr   = new Bio::EnsEMBL::Transcript;
$tr1   = new Bio::EnsEMBL::Transcript;
print "ok 3\n";   
$ex1   = new Bio::EnsEMBL::Exon;
$ex2   = new Bio::EnsEMBL::Exon;
$ex3   = new Bio::EnsEMBL::Exon;
print "ok 4\n";   

$seq = Bio::Seq->new( -id => 'Contig-1' , -seq => 'ATGGCGGATGTTTATGTGGGTGGCCCGGGG' );
$seq2 = Bio::Seq->new( -id => 'Contig-2' , -seq => 'TCAGAAATTTGGGTGTTTTGGCCCTGGTGGTTTGGGTTT' );


$ex1->id("dummy_id_1");
$ex1->contig_id("c_id_1");
$ex1->phase(0);
$ex1->start(8);
$ex1->end(13);
$ex1->attach_seq($seq);
$ex1->strand(1);

$ex2->id("dummy_id_2");
$ex2->contig_id("c_id_2");
$ex2->phase(0);
$ex2->start(18);
$ex2->end(23);
$ex2->attach_seq($seq2);
$ex2->strand(1);

$ex3->id("dummy_id_3");
$ex3->contig_id("c_id_2");
$ex3->phase(0);
$ex3->start(26);
$ex3->end(28);
$ex3->attach_seq($seq2);
$ex3->strand(1);



$tr->add_Exon($ex1);
$tr->add_Exon($ex2);
$trans = Bio::EnsEMBL::Translation->new();
$trans->start_exon_id('dummy_id_1');
$trans->start(1);
$trans->end_exon_id('dummy_id_2');
$trans->end(6);
$tr->translation($trans);

$tr1->add_Exon($ex1);
$tr1->add_Exon($ex2);
$tr1->add_Exon($ex3);
$trans = Bio::EnsEMBL::Translation->new();
$trans->start_exon_id('dummy_id_1');
$trans->start(1);
$trans->end_exon_id('dummy_id_3');
$trans->end(2);
$tr1->translation($trans);

$gene->add_Transcript($tr);
$gene->add_Transcript($tr1);

print "ok 5\n";


$count = 0;
foreach $tr ( $gene->each_Transcript() ) {
	foreach $x ( $tr->each_Exon() ) {
		$count++;
		}
	}

$x =0; # stop it whining...


if( $count != 5 ) {
    print "not ok 6\n";
} else {
    print "ok 6\n";
}


@exons = $gene->each_unique_Exon();

if( scalar @exons != 3 ) {
      print "not ok 7\n";
} else {
       print "ok 7\n";
}

@contigs = $gene->unique_contig_ids();

if( scalar @contigs != 2 ) {
      print "not ok 8\n";
} else {
       print "ok 8\n";
}
	  
foreach $trans ( $gene->each_Transcript ) {
	$pep = $trans->translate();
	}

print "ok 9\n";
$pep = 0;


