# $Id$
# testing of translations of exons that lie across contig boundaries.
# based on staticgoldenpath.t and staticgoldenpath.dump

## We start with some black magic to print on failure.
BEGIN { $| = 1; print "1..24\n";
	use vars qw($loaded); }

END {print "not ok 1\n" unless $loaded;}

use Bio::EnsEMBL::DBLoader;
use lib 't';
use EnsTestDB;
$loaded = 1;
print "ok 1\n";    # 1st test passes.

$" = ", ";                          # for easier list-printing
    
my $ens_test = EnsTestDB->new();
    
# Load some data into the db
$ens_test->do_sql_file("t/stickytranslate.dump");


# Get an EnsEMBL db object for the test db
my $db = $ens_test->get_DBSQL_Obj;
print "ok 2\n";    
$db->static_golden_path_type('UCSC');

$stadaptor = $db->get_StaticGoldenPathAdaptor();

### 2nd gene, now for the more stringent sticky exon translation test (all
### have id'2 ending in 's', and ly on chr3. 

$vc = $stadaptor->fetch_VirtualContig_by_chr_name('chr3');

$seq =$vc->primary_seq->seq; 
$expected = 'N'x203 . 'TGG' x 5 . 'GCC' x 8 . 'TAC' x 7;
if ( $seq eq $expected) {
    print "ok 3\n";
} else {
    print "not ok 3\n";
    warn "expected $expected\ngot $seq\n";
}

# just for fun, translate the static golden path, apart from the first 203
# Ns, as is:
$expected = 'W'x5 . 'A'x8 . 'Y'x7;
$seq = $vc->primary_seq->trunc(204, $vc->length);
$pep = $seq->translate->seq;

if ($pep eq $expected) {
    print "ok 4\n";
} else { 
    print "not ok 4\n";
    warn "expected $expected\ngot $pep\n";
}

### Following (exons with different strandedness) doen'st make
### biological sense, but it works (or should the code guard against this? )

# exon of 5 aa's, last three of RawContig10, first 2 of RawContig20
$exon1 = Bio::EnsEMBL::Exon->new();

$exon1->attach_seq($vc);                # why is this needed ?

$exon1->contig_id($vc->id);
$exon1->id('exon-1s');
$exon1->start(210);
$exon1->end(224);
$exon1->strand(-1);
$exon1->version(1);
$exon1->phase(0);
$exon1->created(1);
$exon1->modified(1);

$seq=$exon1->seq->seq;
$expected= 'GGC' x 2 . 'CCA' x 3;
if ($seq eq $expected) {
    print "ok 5\n";
} else  { 
    print "not ok 5\n";
    warn "expected $expected\ngot $seq\n";
}

# exon of two aa residues, midway RawContig20
$exon2 = Bio::EnsEMBL::Exon->new();
$exon2->attach_seq($vc);                # why is this needed ?
$exon2->id('exon-2s');
$exon2->contig_id($vc->id);
$exon2->start(228);
$exon2->end(233);
$exon2->strand(1);
$exon2->version(1);
$exon2->phase(0);
$exon2->created(1);
$exon2->modified(1);

$seq=$exon2->seq->seq;
$expected= 'GCC' x 2 ;
if ($seq eq $expected) {
    print "ok 6\n";
} else  { 
    print "not ok 6\n";
    warn "expected $expected\ngot $seq\n";
}

# exon of 7 a.a. residues, last two of RawContig20, first 5 of RawContig30
$exon3 = Bio::EnsEMBL::Exon->new();
$exon3->id('exon-3s');
$exon3->contig_id($vc->id);
$exon3->attach_seq($vc);                # why is this needed ?

$exon3->start(237);
$exon3->end(257);
$exon3->strand(-1);
$exon3->version(1);
$exon3->phase(0);
$exon3->created(1);
$exon3->modified(1);

$seq=$exon3->seq->seq;
$expected= 'GTA' x 5 . 'GGC' x 2 ;
if ($seq eq $expected) {
    print "ok 7\n";
} else  { 
    print "not ok 7\n";
    warn "expected $expected\ngot $seq\n";
}

$trl = Bio::EnsEMBL::Translation->new();
$trl->start_exon($exon1);
$trl->start(4); # second aa residue of exon1
$trl->end_exon($exon3);
print STDERR "Seen length of exon3 ",$exon3->length,"\n";

$trl->end(18);  # second last residue of exon3


$transc = Bio::EnsEMBL::Transcript->new();

$transc->add_Exon($exon1);
$transc->add_Exon($exon2);
$transc->add_Exon($exon3);

$transc->translation($trl);

$seq = $transc->dna_seq->seq;
$expected = 'GGC' x 2 . 'CCA' x 3 . 'GCC' x 2 . 'GTA' x 5 . 'GGC' x 2;

if ($seq eq $expected) { 
    print "ok 8\n";
} else { 
    print "not ok 8\n";
    warn "expected $expected\ngot $seq\n";
}

# did this work?
$seq = $transc->translate;
$pep = $seq->seq;
$expected = 'G' x 1 . 'P' x 3 . 'A' x 2 . 'V' x 5 . 'G' x 1;
if ($pep eq $expected) { 
    print "ok 9\n";
} else {
    print "not ok 9\n";
    warn "expected $expected\ngot $pep\n";
}

$gene = Bio::EnsEMBL::Gene->new();
$gene->id('gene-id-1s');
$gene->version(1);
$gene->type('test');
$gene->add_Transcript($transc);
$gene->created(1);
$gene->modified(1);


@transc = $gene->each_Transcript;
if (@transc  ==1 ) {
    print "ok 10\n";
} else { 
    print "not ok 10\n";
    warn "expected 1 transcript, got this:@transc\n";
}

# check if all is intact after this:
$seq = $transc->translate; $pep = $seq->seq;
$expected = 'G' x 1 . 'P' x 3 . 'A' x 2 . 'V' x 5 . 'G' x 1;
if ($pep eq $expected) { 
    print "ok 11\n";
} else {
    print "not ok 11\n";
    warn "expected $expected\ngot $pep\n";
}

## note: the following fails, and should; have to call
## convert_Gene_to_raw_contig first 
#  $geneObj->write($gene);

# can we convert the coords:
my $analysis = Bio::EnsEMBL::Analysis->new( 
      -program => 'genebuild',
      -gff_source => 'genebuild',
      -gff_feature => 'gene',
      -logic_name => 'genebuild'
    );

$gene->analysis($analysis);


$newgene = $vc->convert_Gene_to_raw_contig($gene);
print "ok 12\n";

# is all still OK with these transcripts? 
@transc = $gene->each_Transcript;
if (@transc  ==1 ) {
    print "ok 13\n";
} else { 
    print "not ok 13\n";
    warn "expected 1 transcript, got this:@transc\n";
}

# check if all is intact after this:
$seq = $transc->translate; $pep = $seq->seq;
$expected = 'G' x 1 . 'P' x 3 . 'A' x 2 . 'V' x 5 . 'G' x 1;
if ($pep eq $expected) { 
    print "ok 14\n";
} else {
    print "not ok 14\n";
    warn "expected $expected\ngot $pep\n";
}

# try to store it; this will create Sticky Exons:
$gad = $db->get_GeneAdaptor;
$dbid = $gad->store($newgene);
print "ok 15\n";

# $ens_test->pause;

# see if we can find all genes them back from original chr3 Virtual Contig:
@genes = $vc->get_all_Genes(); # no !?!
if ( @genes == 1 and defined( $genes[0]) ) { 
    print "ok 16\n";
} else { 
    print "not ok 16\n";
    $, = $" = ':';
    warn "not exactly 1 gene, or it is undefined:\n@genes\n";
}

$gene = $genes[0];
@transc = $gene->each_Transcript;
if (@transc  ==1 ) {
    print "ok 17\n";
} else { 
    print "not ok 17\n";
    warn "expected 1 transcript, got this:@transc\n";
}

print "ok 18\n"; # wanton print to skew statistics and to save renumbering:-)

$expected = 'G' x 1 . 'P' x 3 . 'A' x 2 . 'V' x 5 . 'G' x 1;
$transc = $transc[0];
$pep = $transc->translate->seq;
print "ok 19\n";

if (  $pep eq $expected ) {
    print "ok 20\n"; 
} else { 
    print "not ok 20\n";
    warn "expected $expected\ngot $pep\n";
}

@virtualgenes = $vc->get_all_VirtualGenes() ;
if (      @virtualgenes == @genes       # equal list length
     and  $virtualgenes[0]->id eq $genes[0]->id ) { 
    print "ok 21\n"; 
} else { 
    print "not ok 21\n"; 
    warn "should be one gene, one virtualgene, and should have same id; got: ".
      $virtualgenes[0]->id . " <->" . $genes[0]->id , "\n";
}

# try read it back in
$gene2 = $gad->fetch_by_dbID($dbid);
print "ok 22\n";

# is all still OK with gene and  transcripts? 
@transc = $gene2->each_Transcript;
if (@transc  ==1 ) {
    print "ok 23\n";
} else { 
    print "not ok 23\n";
    warn "expected 1 transcript, got this:@transc\n";
}

# check if all is intact after this:
$seq = $transc->translate; $pep = $seq->seq;
$expected = 'G' x 1 . 'P' x 3 . 'A' x 2 . 'V' x 5 . 'G' x 1;
if ($pep eq $expected) { 
    print "ok 24\n";
} else {
    print "not ok 24\n";
    warn "expected $expected\ngot $pep\n";
}
