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
BEGIN { $| = 1; print "1..18\n"; 
	use vars qw($loaded); }

END {print "not ok 1\n" unless $loaded;}

use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DBLoader;
use lib 't';
use EnsTestDB;
$loaded = 1;
print "ok 1\n";    # 1st test passes.
    
my $ens_test = EnsTestDB->new();
    
# Load some data into the db
$ens_test->do_sql_file("t/staticgoldenpath.dump");

# Get an EnsEMBL db object for the test db
my $db = $ens_test->get_DBSQL_Obj;
print "ok 2\n";    
$db->static_golden_path_type('UCSC');

$stadaptor = $db->get_StaticGoldenPathAdaptor();

@array = $stadaptor->fetch_RawContigs_by_fpc_name('ctg123');
if( scalar(@array) != 3 ) {
   print "not ok 3\n";
} else {
   print "ok 3\n";
}

$rc1 = shift @array;
if( $rc1->id ne 'contig1' ) {
   print "not ok 4\n";
} else {
   print "ok 4\n";
}

@array = $stadaptor->fetch_RawContigs_by_chr_name('chr2');
if( scalar(@array) != 3 ) {
   print "not ok 5\n";
} else {
   print "ok 5\n";
}

$rc1 = shift @array;
if( $rc1->id ne 'contig1' ) {
   print "not ok 6\n";
} else {
   print "ok 6\n";
}


$vc = $stadaptor->fetch_VirtualContig_by_fpc_name('ctg123');



#my $vseq=$vc->primary_seq->subseq(10,20);
#print STDERR "Subseq from 2 to 3: ".$vseq->subseq(2,3)."\n";
#if( $vc->primary_seq->seq eq 'AAATTT' ) {
#    print "ok 7\n";
#} else {
#    print "not ok 7\n";
#}
print "ok 7\n";


$vc2 = $stadaptor->fetch_VirtualContig_by_chr_name('chr2');

print "ok 8\n";

#$vc2->_dump_map(\*STDERR);

# ok. lets test some converts

($start,$end,$strand) = $vc2->_convert_start_end_strand_vc('contig1',5,7,1);

if( $start != 203 || $end != 205 || $strand != 1 ) {
    print "not ok 9\n";
    print STDERR "Got $start:$end:$strand\n";
} else {
    print "ok 9\n";
}

($start,$end,$strand) = $vc2->_convert_start_end_strand_vc('contig2',5,7,1);

if( $start != 368 || $end != 370 || $strand != -1 ) {
    print "not ok 10\n";
    print STDERR "Got $start:$end:$strand\n";
} else {
    print "ok 10\n";
}

($rc,$rc_pos,$rc_strand) = $vc2->_vmap->raw_contig_position(203,1);

if( $rc->id ne 'contig1' || $rc_pos != 5 || $rc_strand != 1 ) {
    print "not ok 11\n";
} else {
    print "ok 11\n";
}

($rc,$rc_pos,$rc_strand) = $vc2->_vmap->raw_contig_position(205,-1);

if( $rc->id ne 'contig1' || $rc_pos != 7 || $rc_strand != -1 ) {
    print "not ok 12\n";
} else {
    print "ok 12\n";
}

($rc,$rc_pos,$rc_strand) = $vc2->_vmap->raw_contig_position(368,1);

if( $rc->id ne 'contig2' || $rc_pos != 7 || $rc_strand != -1 ) {
    print "not ok 13\n";
} else {
    print "ok 13\n";
}

($rc,$rc_pos,$rc_strand) = $vc2->_vmap->raw_contig_position(370,-1);

if( $rc->id ne 'contig2' || $rc_pos != 5 || $rc_strand != 1 ) {
    print "not ok 14\n";
} else {
    print "ok 14\n";
}

### write gene test


$gene = Bio::EnsEMBL::Gene->new();
$gene->id('gene-id-1');
$gene->version(1);
$trans = Bio::EnsEMBL::Transcript->new();
$trans->id('trans-id-1');
$trans->version(1);

$dbl = Bio::Annotation::DBLink->new();
$dbl->database('embl');
$dbl->primary_id('AC000012');
$trans->add_DBLink($dbl);

$dbl = Bio::Annotation::DBLink->new();
$dbl->database('swissprot');
$dbl->primary_id('P000012');
$gene->add_DBLink($dbl);

$trl = Bio::EnsEMBL::Translation->new();
$trl->start_exon_id('exon-1');
$trl->end_exon_id('exon-2');
$trl->start(203);
$trl->end(370);
$trl->id('trl-id');
$trl->version(1);
$trans->translation($trl);
$trans->created(1);
$trans->modified(1);
$gene->add_Transcript($trans);
$gene->created(1);
$gene->modified(1);

$exon = Bio::EnsEMBL::Exon->new();
$exon->id('exon-1');
$exon->start(200);
$exon->end(205);
$exon->strand(1);
$exon->version(1);
$exon->phase(0);
$exon->created(1);
$exon->modified(1);
$exon->contig_id($vc2->id);
$trans->add_Exon($exon);
$exon{'exon-1'} = $exon;


$sf = Bio::EnsEMBL::FeatureFactory->new_feature_pair();

$sf->start(200);
$sf->end(205);
$sf->hstart(100);
$sf->hend(110);
$sf->strand(1);
$sf->seqname($vc2->id);
$sf->hseqname('other');
$sf->score(100);
$sf->primary_tag('similarity');
$sf->source_tag('someone');
$sf->feature2->primary_tag('similarity');
$sf->feature2->source_tag('someone');

$analysis = Bio::EnsEMBL::FeatureFactory->new_analysis();
$analysis->program('program');
$analysis->program_version('version-49');
$analysis->gff_source('source');
$analysis->gff_feature('feature');

$sf->analysis($analysis);
$sf->feature2->analysis($analysis);
$sf->hstrand(1);
$sf->hscore(100);

$exon->add_Supporting_Feature($sf);

$exon = Bio::EnsEMBL::Exon->new();
$exon->id('exon-2');
$exon->start(367);
$exon->end(373);
$exon->strand(1);
$exon->version(1);
$exon->phase(0);
$exon->created(1);
$exon->modified(1);
$exon->contig_id($vc2->id);
$trans->add_Exon($exon);
$exon{'exon-2'} = $exon;

$exon = Bio::EnsEMBL::Exon->new();
$exon->id('exon-3');
$exon->start(373);
$exon->end(380);
$exon->strand(1);
$exon->version(1);
$exon->phase(0);
$exon->created(1);
$exon->modified(1);
$exon->contig_id($vc2->id);
$trans->add_Exon($exon);


$newgene = $vc2->convert_Gene_to_raw_contig($gene);
print "ok 15\n";

$db->write_Gene($newgene);

print "ok 16\n";

# retrieve gene, check it is on the correct coordinates

$newgene = $db->gene_Obj->get('gene-id-1');
($trans) = $newgene->each_Transcript;
if( !defined $trans ) { die "Didn't even get a transcript!"; }

($exon1,$exon2,$exon3) = $trans->each_Exon;

if( $exon1->id ne 'exon-1' || $exon1->start != 2 ||
    $exon1->end != 7 || $exon1->strand != 1 ||
    $exon2->id ne 'exon-2' || $exon2->start != 2 ||
    $exon2->end != 8 || $exon2->strand != -1 ) {
    print "not ok 17\n";
    print STDERR "exon1 ",$exon1->start,":",$exon1->end,";",$exon1->strand," ",
    "exon2 ",$exon2->start,":",$exon2->end,";",$exon2->strand,"\n";
} else {
    print "ok 17\n";
}

# exon3 is a sticky exon, starting 1 ending 8 on a sticky-wicket
if( $exon3->start != 1 || $exon3->end != 8 || $exon3->seq->seq ne 'GGTTTTTT' ) {
    print "not ok 18\n";
    print STDERR "Exon3 ",$exon3->seqname,"  ",$exon3->start,":",$exon3->end,";",$exon3->strand," ",$exon3->seq->seq,"\n";
} else {
    print "ok 18\n";
}







