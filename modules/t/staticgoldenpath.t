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
BEGIN { $| = 1; print "1..33\n"; 
	use vars qw($loaded); }

END {print "not ok 1\n" unless $loaded;}

use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DBLoader;
#use EnsemblExt;

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


#$vc->_dump_map(\*STDERR);


my $vseq=$vc->primary_seq;
print "ok 7\n";

# for more serious vseq tests, read below...

$vc2 = $stadaptor->fetch_VirtualContig_by_chr_name('chr2');

print "ok 8\n";

#$vc2->_dump_map(\*STDERR);

# ok. lets test some converts

($start,$end,$strand) = $vc2->_convert_start_end_strand_vc('contig1',5,7,1);

if( $start != 204 || $end != 206 || $strand != 1 ) {
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

($rc,$rc_pos,$rc_strand) = $vc2->_vmap->raw_contig_position(204,1);

if( $rc->id ne 'contig1' || $rc_pos != 5 || $rc_strand != 1 ) {
    print "not ok 11\n";
} else {
    print "ok 11\n";
}

($rc,$rc_pos,$rc_strand) = $vc2->_vmap->raw_contig_position(206,-1);

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
$gene->type('test');
$trans = Bio::EnsEMBL::Transcript->new();
$trans->id('trans-id-1');
$trans->version(1);

my $analysis = Bio::EnsEMBL::Analysis->new( 
      -program => 'genebuild',
      -gff_source => 'genebuild',
      -gff_feature => 'gene',
      -logic_name => 'genebuild'
    );

$gene->analysis($analysis);

$exDB = Bio::EnsEMBL::DBEntry->new
    ( -primary_id => 'AC000012',
      -display_id => 'ACC_optional',
      -version => 1,
      -release => 2,
      -dbname => 'embl' );


$trans->add_DBLink($exDB);

#$dbl = Bio::Annotation::DBLink->new();
#$dbl->database('swissprot');
#$dbl->primary_id('P000012');

$exDB = Bio::EnsEMBL::DBEntry->new
    ( -primary_id => 'P000012',
      -display_id => 'ACC_optional',
      -version => 1,
      -release => 2,
      -dbname => 'swissprot' );


$gene->add_DBLink($exDB);

$gene->add_Transcript($trans);
$gene->created(1);
$gene->modified(1);

$exon = Bio::EnsEMBL::Exon->new();
$exon->id('exon-1');
$exon->start(201);
$exon->end(206);
$exon->strand(1);
$exon->version(1);
$exon->phase(0);
$exon->created(1);
$exon->modified(1);
$exon->contig_id($vc2->id);
$trans->add_Exon($exon);


$sf = Bio::EnsEMBL::FeatureFactory->new_feature_pair();

$sf->start(201);
$sf->end(206);
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

my $analysis = Bio::EnsEMBL::Analysis->new( 
      -program => 'supporting',
      -gff_source => 'sf',
      -gff_feature => 'yadda',
      -logic_name => 'genewise'
    );

$sf->analysis($analysis);
$sf->feature2->analysis($analysis);
$sf->hstrand(1);
$sf->hscore(100);

$exon->add_Supporting_Feature($sf);

$trl = Bio::EnsEMBL::Translation->new();
$trl->start_exon($exon);
$trl->start(3);
$trl->end(4);
$trl->id('trl-id');
$trl->version(1);

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
$trl->end_exon($exon);

$trans->translation($trl);
$trans->created(1);
$trans->modified(1);

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
for $ex ( $newgene->get_all_Exons() ) {
	print STDERR "Exon ",$exon->contig_id,"\n";
}

$db->write_Gene($newgene);

#$ens_test->pause;

print "ok 16\n";

# retrieve gene, check it is on the correct coordinates

$newgene = $db->get_GeneAdaptor->fetch_by_dbID($newgene->dbID);


($trans) = $newgene->each_Transcript;
if( !defined $trans ) { die "Didn't even get a transcript!"; }

($exon1,$exon2,$exon3) = $trans->each_Exon;

if( $exon1->start != 2 ||
    $exon1->end != 7 || $exon1->strand != 1 ||
    $exon2->start != 2 ||
    $exon2->end != 8 || $exon2->strand != -1 ) {
    print "not ok 17\n";
    print STDERR "exon1 ",$exon1->start,":",$exon1->end,";",$exon1->strand," ",
    "exon2 ",$exon2->start,":",$exon2->end,";",$exon2->strand,"\n";
} else {
    print "ok 17\n";
}

# exon3 is a sticky exon, starting 1 ending 8 on a sticky-wicket
print STDERR "Exon 3 has ",$exon3->start," ",$exon3->end,"\n";

if( $exon3->start != 1 || $exon3->end != 8 || $exon3->seq->seq ne 'GGTTTTTT' ) {
    print "not ok 18\n";
    print STDERR "Exon3 ",$exon3->seqname,"  ",$exon3->start,":",$exon3->end,";",$exon3->strand," ",$exon3->seq->seq,"\n";
} else {
    print "ok 18\n";
}


if( $vseq->subseq(1,10) ne 'AAAAAAAAAT' ) {
    print "not ok 19\n";
    print STDERR "Sequence ",$vseq->subseq(1,10),"\n";
} else {
    print "ok 19\n";
}

if( $vseq->subseq(1,40) ne 'AAAAAAAAATTTTTTTTTTAAAAAAAAAATTTTTTNNNNN' ) {
    print "not ok 20\n";
    print STDERR "Sequence ",$vseq->subseq(1,40),"\n";
} else {
    print "ok 20\n";
}

#print STDERR "Whole is ",$vseq->seq,"\n";


# testing inversions

my $invert = $vc2->invert();

if( $invert->length != $vc2->length ) {
    print "not ok 21\n";
} else {
    print "ok 21\n";
}


if( $invert->seq ne $vc2->revcom->seq ) {
    print "not ok 22\n";
} else {
    print "ok 22\n";
}

my @fp = $invert->get_all_SimilarityFeatures();
my $fp = shift @fp;

if( $fp->start != 9 || $fp->end != 19 ) {
    print "not ok 23\n";
    print STDERR "Feature ",$fp->start," ",$fp->end,"\n";
} else {
    print "ok 23\n";
}

#$vc2->_dump_map(\*STDERR);
#print STDERR "//\n";
#$invert->_dump_map(\*STDERR);

# testing get on vc is good
$vcsave = $vc2;
@genes = $vc2->get_all_Genes();
$gene = shift @genes;
if( !defined $gene ) {
    print "not ok 24\n";
} else {
    $error = 0;
    foreach $exon ( $gene->each_unique_Exon ) {
	if( $exon->seqname ne $vc2->id ) {
	    print STDERR "Got exon on ",$exon->seqname," not ",$vc2->id,"\n";
	    $error = 1;
	}
	#print STDERR "Got exon [$exon] on ",$exon->seqname," not ",$vc2->id,"\n";
    }
    if( $error == 0 ) {
	print "ok 24\n";
    } else {
	print "not ok 24\n";
    }
}


my $vc = $stadaptor->fetch_VirtualContig_by_chr_start_end('chr2',1,380);
print "ok 25\n";

if( $vc->length eq length($vc->seq) && $vc->length == 380) {
    print "ok 26\n";
} else {
    print "not ok 26\n";
}



$vc = $stadaptor->fetch_VirtualContig_by_chr_start_end('chr2',250,350);


print "ok 27\n";

if( $vc->length eq length($vc->seq) && $vc->length == 101) {
    print "ok 28\n";
} else {
    print "not ok 28\n";
}

$vc = $stadaptor->fetch_VirtualContig_by_chr_start_end('chr2',250,260);

print "ok 29\n";

if( $vc->length eq length($vc->seq) && $vc->length == 11) {
    print "ok 30\n";
} else {
    print "not ok 30\n";
}


# my $vc = $stadaptor->fetch_VirtualContig_by_clone('pog',200);

#print "ok 31\n";


#my $vc = $stadaptor->fetch_VirtualContig_by_contig('contig2',200);
#print "ok 32\n";

$contig = $stadaptor->fetch_VirtualContig_by_chr_start_end('chr2',10,300);

if ($contig->isa(Bio::EnsEMBL::Virtual::StaticContig)){print "ok 31\n";}
else {print "not ok 31\n";}


@genes = $vcsave->get_all_Genes_exononly();

#foreach $g ( @genes ) {
#	print STDERR "Got ",$g->id,"\n";
#	foreach $e ( $g->each_unique_Exon ) {
#		print STDERR " Got exon ",$e->id,"\n";
#	}
#}

if( scalar(@genes) == 2 ) {
	print "ok 32\n";
}

# testing C-extensions

#$vcsave->_use_cext_get(1);
#@f = $vcsave->get_all_SimilarityFeatures_above_score('swir',80);
#
#foreach $f ( @f ) {
#   print "Got $f ",$f->start," ",$f->end,"\n";
#}


$chadp = $db->get_ChromosomeAdaptor();

$chadp->fetch_by_chrname('chr2');

print "ok 33\n";

#my @m = $chr->get_landmark_MarkerFeatures();
#
#if( scalar(@m) > 0 ) {
#   print "ok 34\n";
#}

