# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

use strict;
use warnings;

use Test::More;
use Test::Warnings qw(warning);

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Intron;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

ok( $multi );

my $db = $multi->get_DBAdaptor('core' );
my $ontology = Bio::EnsEMBL::Test::MultiTestDB->new('ontology');
my $odb;
warning { $odb = $ontology->get_DBAdaptor("ontology"); };

my $sa = $db->get_SliceAdaptor();

my $slice = $sa->fetch_by_region('chromosome', "20", 30_249_935, 31_254_640 );

my $genes = $slice->get_all_Genes();

my $translates = 1;
my $utr_trans;

note( "Checking if all Transcripts translate and transform" );

for my $gene ( @$genes ) {
  for my $trans ( @{$gene->get_all_Transcripts()} ) {
    if( $trans->translate()->seq() =~ /\*./ ) {
      $translates = 0;
      note( $trans->stable_id()." does not translate." );
      last;
    }
    if( $trans->coding_region_start() != $trans->start() &&
	$trans->coding_region_end() != $trans->end() ) {
      $utr_trans = $trans->stable_id();
    }
  }

  $gene = $gene->transform( "contig" );
  next if( ! $gene );	
  for my $trans ( @{$gene->get_all_Transcripts()} ) {
    if( $trans->translate()->seq() =~ /\*./ ) {
      $translates = 0;
      note( $trans->stable_id()." does not translate." );
      last;
    }
  }
}

note( "utr Transcript is $utr_trans" );

ok( $translates );

#
# Try pulling off genes from an NTContig and making sure they still translate.
# This is a fairly good test of the chained coordinate mapping since the
# transcripts are stored in chromosomal coordinates and there is no direct
# mapping path to the supercontig coordinate system.
#
my $supercontig = $sa->fetch_by_region('supercontig', "NT_028392");

my $transcripts = $supercontig->get_all_Transcripts();

note( "Checking if all transcripts on NTContig NT_028392 translate and transform" );

$translates = 1;
for my $trans ( @$transcripts ) {
  if( $trans->translate()->seq() =~ /\*./ ) {
    $translates = 0;
    note( $trans->stable_id()." does not translate." );
    last;
  }
  note($trans->stable_id() .  ":" . $trans->translate()->seq() . "\n");
}

ok($translates);


my $ta = $db->get_TranscriptAdaptor();

note ("Transcript->list_dbIDs");
my $ids = $ta->list_dbIDs();
ok (@{$ids});

note ("Transcript->list_stable_ids");
my $stable_ids = $ta->list_stable_ids();
ok (@{$stable_ids});

my $tr = $ta->fetch_by_display_label( "DNMT3B" );
is($tr->dbID, 21737, 'Fetched correct dbID');


$tr = $ta->fetch_by_stable_id( "ENST00000217347" );

$tr = $tr->transform('contig');

ok( $tr );

note ( "External transcript name: " . $tr->external_name );
is ($tr->external_name, "MAPRE1", 'Fetched correct external name');

note ( "External transcript dbname: " . $tr->external_db );
is ( $tr->external_db, 'HUGO' , 'Fetched correct external db');

note ( "Display_xref_id: " . $tr->display_xref->dbID() );
is ( $tr->display_xref->dbID(), 97759, 'Fetched correct display xref id' );
ok( test_getter_setter( $tr, "display_xref", 42 ));

ok( test_getter_setter( $tr, "dbID", 100000 ));
ok( test_getter_setter( $tr, "biotype", "NOVEL" ));
ok( test_getter_setter( $tr, "created_date", time() ));
ok( test_getter_setter( $tr, "modified_date", time() ));


my @date_time = localtime( $tr->created_date());
ok( $date_time[3] == 6 && $date_time[4] == 11 && $date_time[5] == 104 );

@date_time = localtime( $tr->modified_date());
ok( $date_time[3] == 6 && $date_time[4] == 11 && $date_time[5] == 104 );


ok( $tr->translation->isa( "Bio::EnsEMBL::Translation" ));

is( $tr->start(), 79874, 'Start is correct' );

is( $tr->end(), 110306, 'End is correct' );

is ( substr( $tr->spliced_seq(), 0, 10 ), "ACGAGACGAA", 'Start of spliced seq is correct' ); 
is ( substr( $tr->spliced_seq(1), 0, 10 ), "acgagacgaa", 'Spliced seq with utr lower casing is correct');
is ( length($tr->spliced_seq()), length($tr->spliced_seq(1)), "Spliced seq with or without utr lower casing has the same length");
is ( $tr->spliced_seq(), uc($tr->spliced_seq(1)), "Spliced seq is identical to upper case utr masked spliced seq");
is ( substr($tr->spliced_seq(1), 61, 6), 'aagATG', 'Start mask boundary on forward stand transcript is correct' );
is ( substr($tr->spliced_seq(1), 865, 6), 'TATtaa', 'End mask boundary on forward stand transcript is correct' );

is ( substr( $tr->translateable_seq(),0,10 ), "ATGGCAGTGA", 'Start of translateable sequence is correct' );

is( $tr->coding_region_start(), 85834, 'Correct coding region start' );

is( $tr->coding_region_end(), 108631, 'Correct coding region end' );

my @pepcoords = $tr->pep2genomic( 10, 20 );
is( $pepcoords[0]->start(), 85861, 'Correct translation start' );

my $t_start = $tr->start;
my $t_end   = $tr->end;
my $t_strand = $tr->get_all_Exons->[0]->strand;
my $pep_len = ($tr->cdna_coding_end - $tr->cdna_coding_start + 1) / 3;

my $coord_num = 0;
my $coding_region_start = $tr->coding_region_start;
my $coding_region_end   = $tr->coding_region_end;
#expect coordinate for each exon with coding sequence 
foreach my $e (@{$tr->get_all_Exons}) {
  if($e->end > $coding_region_start && $e->start < $coding_region_end) {
    $coord_num++;
  }
}

my @coords = ();
my @gaps = ();
my $coord_len = 0;
foreach my $c ($tr->genomic2pep($t_start, $t_end, $t_strand)) {
  if($c->isa('Bio::EnsEMBL::Mapper::Gap')) {
    push @gaps, $c;
  } else {
    push @coords, $c;
  }
}

is(scalar(@coords), $coord_num, 'Number of coords is correct');
my ($last_coord) = reverse(@coords); 
is($pep_len, $last_coord->end, 'Peptide length matched end coordinate');

note( "start Exon: ".$tr->start_Exon->stable_id() );
note( "end Exon: ".$tr->end_Exon->stable_id() );

is($tr->cdna_coding_start, 65, 'Correct cdna coding start');
ok(test_getter_setter($tr, 'cdna_coding_start', 99));

is( substr( $tr->five_prime_utr()->seq(), -5, 5), "CGAAG", 'Five prime utr seq is correct' ); 

is($tr->cdna_coding_end, 868, 'Correct cdna coding end');
ok(test_getter_setter($tr, 'cdna_coding_end', 102));

is( substr( $tr->three_prime_utr()->seq(), -5, 5  ), "TTCAA", 'Three prime utr seq is correct');

is( scalar( @{$tr->get_all_Exons()} ), 7, 'Transcript has 7 exons' );

$tr->flush_Exons();

is( scalar( @{$tr->get_all_Exons()} ), 0, 'No exons left after flushing' );

# Fetch a fresh tr, check incomplete codon behavior
$tr = $ta->fetch_by_stable_id( "ENST00000300425" );

# By default the incomplete codon should be dropped
is( $tr->translate()->seq() =~ /P$/, 1, "Incomplete codon is not translated");
is( $tr->translate(1)->seq() =~ /PL$/, 1, "Incomplete codon is padded then translated");

# get a fresh tr to check the update method
$tr = $ta->fetch_by_stable_id( "ENST00000217347" );

$multi->save('core', 'transcript', 'meta_coord');

# the first update should have no effect
$ta->update($tr);

my $up_tr = $ta->fetch_by_stable_id( "ENST00000217347" );
is( $up_tr->display_xref->dbID(), 97759, 'Fetched the correct dbID' );

my $dbentryAdaptor = $db->get_DBEntryAdaptor();

$tr->display_xref($dbentryAdaptor->fetch_by_dbID( 614 ));
$ta->update($tr);

$up_tr = $ta->fetch_by_stable_id( "ENST00000217347" );
is ( $up_tr->display_xref->dbID(), 614, 'Fetched the correct display xref id');

$multi->restore('core', 'transcript', 'meta_coord');

#
# Test spliced_seq on a reverse strand transcript
#

$tr = $ta->fetch_by_stable_id( "ENST00000246229" );

is ( substr( $tr->spliced_seq(), 0, 10 ), "ATGGCCCGAC", 'Start of spliced seq is correct, rev strand' );
is ( substr( $tr->spliced_seq(1), 0, 10 ), "atggcccgac", 'Spliced seq with utr lower casing is correct, rev strand');
is ( length($tr->spliced_seq()), length($tr->spliced_seq(1)), "Spliced seq with or without utr lower casing has the same length, rev strand");
is ( $tr->spliced_seq(), uc($tr->spliced_seq(1)), "Spliced seq is identical to upper case utr masked spliced seq, rev strand");
is ( substr($tr->spliced_seq(1), 199, 6), 'gccATG', 'Start mask boundary on forward stand transcript is correct, rev strand' );
is ( substr($tr->spliced_seq(1), 1687, 6), 'CAGtag', 'End mask boundary on forward stand transcript is correct, rev strand' );


my $interpro = $ta->get_Interpro_by_transid("ENST00000252021");
foreach my $i (@$interpro) {
  note($i);
}
###currently no interpro info in the test db
ok(@$interpro == 1);

#
# test fetch_all_by_external_name
#

($tr) = @{$ta->fetch_all_by_external_name('PLAGL2')};
is($tr->stable_id, 'ENST00000246229', 'Fetched correct transcript by external name');

#
# test fetch_by_translation_id
#

$tr = $ta->fetch_by_translation_id(21734);

is($tr->stable_id, 'ENST00000201961', 'Fetched correct transcript by translation id');

is($tr->display_id(), $tr->stable_id(), 'Transcript stable id and display id are identical');

#
# test TranscriptAdaptor::fetch_all_by_biotype
#
note("Test fetch_all_by_biotype");
my @transcripts = @{$ta->fetch_all_by_biotype('protein_coding')};
is(@transcripts, 27, 'Fetching all protein coding transcript');
my $transcriptCount = $ta->count_all_by_biotype('protein_coding');
is($transcriptCount, 27, 'Counting all protein coding');
@transcripts = @{$ta->fetch_all_by_biotype(['protein_coding','pseudogene'])};
is(@transcripts, 27, 'Got 27 transcript');
$transcriptCount = $ta->count_all_by_biotype(['protein_coding', 'pseudogene']);
is($transcriptCount, 27, 'Count by biotype is correct');

#
# test TranscriptAdaptor::fetch_all_by_Slice
#
note("Test fetch_all_by_Slice");
@transcripts = @{$ta->fetch_all_by_Slice($slice)};
$transcriptCount = $ta->count_all_by_Slice($slice);
is(@transcripts, $transcriptCount, "Counted as many transcripts as were fetched from slice");

#
# test TranscriptAdaptor::fetch_all_by_source
#
note("Test fetch_all_by_source");
@transcripts = @{$ta->fetch_all_by_source('ensembl')};
note "Got ".scalar(@transcripts)." ensembl transcripts\n";
is(24, scalar(@transcripts));
$transcriptCount = $ta->count_all_by_source('ensembl');
is(24, $transcriptCount);
@transcripts = @{$ta->fetch_all_by_source(['havana','vega'])};
note "Got ".scalar(@transcripts)." (havana, vega) transcripts\n";
is(3, scalar(@transcripts));
$transcriptCount = $ta->count_all_by_source(['havana', 'vega']);
is(3, $transcriptCount);

#
# test TranscriptAdaptor::fetch_all
#
note("Test fetch_all");
@transcripts = @{ $ta->fetch_all() };
is(27, scalar(@transcripts), "Got 27 transcripts");

#
# test TranscriptAdaptor::fetch_all_by_GOTerm
#
note("Test fetch_all_by_GOTerm");
my $go_adaptor;
warning { $go_adaptor = $odb->get_OntologyTermAdaptor(); };
my $go_term = $go_adaptor->fetch_by_accession('GO:0070363');
@transcripts = @{ $ta->fetch_all_by_GOTerm($go_term) };
is(scalar(@transcripts), 0, "Found 0 genes with fetch_all_by_GOTerm");


#
# test TranscriptAdaptor::fetch_all_by_exon_supporting_evidence
#
note("Test fetch_all_by_exon_supporting_evidence");
@transcripts = @{ $ta->fetch_all_by_exon_supporting_evidence('BRCA2', 'dna_align_feature') };
is(scalar(@transcripts), 0, "No transcripts with BRCA2 daf");
@transcripts = @{ $ta->fetch_all_by_exon_supporting_evidence('AK015740.1', 'dna_align_feature') };
is(scalar(@transcripts), 1, "1 transcript with AK015740.1 daf");
@transcripts = @{ $ta->fetch_all_by_exon_supporting_evidence('Q9NUG5', 'protein_align_feature') };
is(scalar(@transcripts), 1, "1 transcript with Q9NUG5 paf");

#
# test TranscriptAdaptor::fetch_all_by_transcript_supporting_evidence
#
note("Test fetch_all_by_transcript_supporting_evidence");
@transcripts = @{ $ta->fetch_all_by_transcript_supporting_evidence('BRCA2', 'dna_align_feature') };
is(scalar(@transcripts), 0, "No transcripts with BRCA2 daf");
@transcripts = @{ $ta->fetch_all_by_transcript_supporting_evidence('AK015740.1', 'dna_align_feature') };
is(scalar(@transcripts), 0, "0 transcripts with AK015740.1 daf");
@transcripts = @{ $ta->fetch_all_by_transcript_supporting_evidence('Q9NUG5', 'protein_align_feature') };
is(scalar(@transcripts), 0, "No transcripts with Q9NUG5 paf");


#
# Test get_all_Introns by joining Exons and introns
# and comparing it to the original
#

foreach my $stable_id (qw(ENST00000201961 ENST00000217347)){ #test both strands
#foreach my $stable_id (qw(ENST00000217347)){ #test both strands

  my $transcript_adaptor = $db->get_TranscriptAdaptor();
  my $transcript = 
    $transcript_adaptor->fetch_by_stable_id($stable_id);


  my @exons = (@{$transcript->get_all_Exons()});
  my @introns = (@{$transcript->get_all_Introns()});

  my @cds = (@{$transcript->get_all_CDS()});
  my @cds_introns = (@{$transcript->get_all_CDS_Introns()});

  my $orig_seq = $transcript->slice->subseq(
					    $transcript->start(),
					    $transcript->end(), 
					    $transcript->strand());

  my $cds_orig_seq = $transcript->slice->subseq(
                                            $transcript->coding_region_start(),
                                            $transcript->coding_region_end(),
                                            $transcript->strand());

  my $idl=0;
  my $new_seq = $exons[0]->seq()->seq();
  foreach my $intron (@introns){
    $new_seq .= $intron->seq;
    $new_seq .= $exons[$idl+1]->seq->seq();
    $idl++;

  }

  is($orig_seq, $new_seq, 'Correct new origin seq');

  my $cds_idl=0;
  my $new_cds_seq = $cds[0]->seq();
  foreach my $cds_intron (@cds_introns){
    $new_cds_seq .= $cds_intron->seq;
    $new_cds_seq .= $cds[$cds_idl+1]->seq();
    $cds_idl++;

  }

  is($cds_orig_seq, $new_cds_seq, 'Correct new cds origin seq');

}


#
# regression test:  five_prime_utr and three_prime_utr were failing
# for transcripts that had no UTR. undef should have been returned instead
#
$tr = $ta->fetch_by_stable_id('ENST00000246203');

my $three_prime = $tr->three_prime_utr();

ok(!defined($three_prime));

my $five_prime = $tr->five_prime_utr();

ok(!defined($five_prime));



#
# test removal of transcript
#
my $tl_count = count_rows($db, "translation");
my $ex_tr_count = count_rows($db, "exon_transcript");
my $tr_count = count_rows($db, "transcript");

my $ex_tr_minus = @{$tr->get_all_Exons()};



$multi->save("core", "transcript", "translation",
             "protein_feature", "exon",
             "exon_transcript", "object_xref",
             "supporting_feature", "dna_align_feature","protein_align_feature",
             'xref', "ontology_xref", "identity_xref", 'meta_coord');

$ta->remove($tr);

ok(!defined($tr->dbID()));
ok(!defined($tr->adaptor()));

is( count_rows( $db, "transcript"), ($tr_count - 1), 'Row count matches transcript count');
is( count_rows( $db, "translation"), ($tl_count - 1), 'Row count matches translation count');
is( count_rows( $db, "exon_transcript"), ($ex_tr_count - $ex_tr_minus), 'Row count matches exon count');

#
# test removal of transcript with supporting changes at the gene level
#

$tr = $ta->fetch_by_stable_id('ENST00000278995');
my $gene = $tr->get_Gene;
print $gene->stable_id."\n\n";
note(join "\n",map { $_->stable_id } @{$gene->get_all_Transcripts});

$ta->remove($tr,1);

ok(! grep { $_->stable_id eq 'ENST00000278995'} @{$gene->get_all_Transcripts});

# note(join "\n",map { $_->stable_id } @{$gene->get_all_Transcripts});
# old coords 30274334  30300924, the old gene above is still pointing to old coordinates, but remove cannot find it to change it.
my $new_copy_gene = $ta->fetch_by_stable_id('ENST00000310998')->get_Gene;
cmp_ok($new_copy_gene->start, '==', 30274334, 'Shortened Gene starts as before');
cmp_ok($new_copy_gene->end, '==', 30298904, 'Shortened Gene ends earlier');

$tr = $ta->fetch_by_stable_id('ENST00000278995');
ok(! defined $tr);
#
# test _rna_edit for transcripts
#

$tr = $ta->fetch_by_stable_id( "ENST00000217347" );

#
# 5 prime UTR editing
#

$tr->edits_enabled(1);

my $seq1 = $tr->spliced_seq();
my $tlseq1 = $tr->translateable_seq();

my $cdna_cds_start1 = $tr->cdna_coding_start();
my $cdna_cds_end1   = $tr->cdna_coding_end();

my $attrib = Bio::EnsEMBL::Attribute->new
  (-code => '_rna_edit',
   -value => "1 6 GATTACA",
   -name => "RNA editing");

$tr->add_Attributes( $attrib );

my $seq2 = $tr->spliced_seq(1);

ok( $seq1 ne $seq2 );
ok( $seq2 =~ /^GATTACA/ );

my $cdna_cds_start2 = $tr->cdna_coding_start();
my $cdna_cds_end2   = $tr->cdna_coding_end();

is($cdna_cds_start1, $cdna_cds_start2 - 1, 'Both cdna starts match');
is($cdna_cds_end1, $cdna_cds_end2   - 1, 'Both cdna ends match');


# insert just at the start of the translation
# makes it longer. (For non phase zero start exons)
# cdna_coding_start for this transcript is 65, 64 is just before that

$attrib = Bio::EnsEMBL::Attribute->new
  ( -code => '_rna_edit',
    -value => "65 64 NNN",
    -name => "RNA editing");

$tr->add_Attributes( $attrib );

my $tlseq2 = $tr->translateable_seq();

ok( $tlseq1 ne $tlseq2 );
ok( $tlseq2 =~ /^NNNATG/ );
ok( $tlseq1 eq substr( $tlseq2,3 ));

is($cdna_cds_start2, $tr->cdna_coding_start(), 'Cds start matches cdna coding start');
is($cdna_cds_end2  , $tr->cdna_coding_end() - 3, 'Coding end -3 matches cds end');

# test that the edits can be disabled
$tr->edits_enabled(0);

is($tr->cdna_coding_start(), $cdna_cds_start1, 'Cds start matches cdna coding start');
is($tr->cdna_coding_end()  , $cdna_cds_end1, 'Cds end matches cdna coding end');

is($tr->spliced_seq(), $seq1, 'Spliced sequence is correct');
is($tr->translateable_seq(), $tlseq1, 'Translateable sequence is correct');

#
# try save and retrieve by lazy load
#

$multi->hide( "core", "transcript_attrib" );
my $attribAdaptor = $db->get_AttributeAdaptor();

$attribAdaptor->store_on_Transcript($tr->dbID, $tr->get_all_Attributes);

$tr = $ta->fetch_by_stable_id( "ENST00000217347" );
$tr->edits_enabled(1);

#print " $tr->translateable_seq() : " .  $tr->translateable_seq() . "\n";
#print "$tr->spliced_seq() : " . $tr->spliced_seq() . "\n";

is( $tr->translateable_seq(), $tlseq2, 'Translateable sequence is correct' );
ok( $tr->spliced_seq() =~ /^GATTACA/ );

$multi->restore();

#
# test that the transcript mapper handles edits
#
$tr = $ta->fetch_by_translation_id(21734);
$slice = $tr->feature_Slice();
$tr = $tr->transfer($slice);
test_trans_mapper_edits($tr);

# test again with reverse strand
$tr = $ta->fetch_by_translation_id(21734);
$slice->invert();
$tr = $tr->transfer($slice);
test_trans_mapper_edits($tr);



# test that mitochondrial transcripts can be corrected translated with
# an alternate codon table

$tr = $ta->fetch_by_stable_id('ENST00000355555');

note($tr->translate->seq());

is($tr->translate->seq(), 'MNFALILMINTLLALLLMIITFWLPQLNGYMEKSTPYECGFDPMSPARVPFSMKFFLVAITFLLFDLEIALLLPLPWALQTTNLPLMVMSSLLLIIILALSLAYEWLQKGLDWAE', 'Translates correctly');



# check that attributes are stored when transcript is stored

$tr = $ta->fetch_by_stable_id( "ENST00000217347" );
my $g = $db->get_GeneAdaptor->fetch_by_transcript_id($tr->dbID());


$tr->translation()->adaptor(undef);
$tr->translation()->dbID(undef);

# unstore the transcript so it can be stored again

foreach my $ex (@{$tr->get_all_Exons()}) {
  $ex->dbID(undef);
  $ex->adaptor(undef);
}

$tr->dbID(undef);
$tr->adaptor(undef);

{
  # testing transform with gaps in introns
  my $tr = $ta->fetch_by_dbID( 21739 );
  my $mapped_tr = $tr->transform( "alt_chrom" );
  is( $tr->spliced_seq(), $mapped_tr->spliced_seq(), 'Transcript seq matches mapped seq' );
}

$multi->hide('core', 'transcript', 'transcript_attrib', 'translation',
             'exon_transcript', 'exon', 'meta_coord');


my $attrib1 = Bio::EnsEMBL::Attribute->new
  ( -code => '_rna_edit',
    -value => "65 64 NNN",
    -name => "RNA editing");

$tr->add_Attributes($attrib1);

my $attrib2 = Bio::EnsEMBL::Attribute->new
  ( -code => '_rna_edit',
    -value => "66 65 NNN",
    -name => "RNA editing");

$tr->add_Attributes($attrib2);

$ta->store($tr, $g->dbID());

is(count_rows($db, 'transcript_attrib'), 2, '2 rows for transcript_attrib');


#
# specifically test genomic2cdna given coords either side of an insertion
#

$tr = $ta->fetch_by_stable_id( "ENST00000217347" );

my $tmp_coding_start = $tr->coding_region_start;

$tr->add_Attributes(
  Bio::EnsEMBL::Attribute->new(
    -code => '_rna_edit',
    -value => "68 67 NNN",
    -name => "RNA editing"
  )
);

$tr->edits_enabled(1);

is(substr($tr->translateable_seq, 0, 6), 'ATGNNN', 'insertion mapper test - edit applied');

is_deeply(
  [$tr->get_TranscriptMapper->genomic2cdna($tmp_coding_start, $tmp_coding_start + 4, 1)],
  [
    bless( {
      'coord_system' => undef,
      'strand' => 1,
      'id' => 'cdna',
      'rank' => 0,
      'end' => 67,
      'start' => 65
    }, 'Bio::EnsEMBL::Mapper::Coordinate' ),
    bless( {
      'coord_system' => undef,
      'strand' => 1,
      'id' => 'cdna',
      'rank' => 0,
      'end' => 72,
      'start' => 71
    }, 'Bio::EnsEMBL::Mapper::Coordinate' )
  ],
  'insertion mapper test - coords either side of inserted block'
);


$multi->restore('core');

#
# tests for translation start/end within a seq_edit
#
$tr = $ta->fetch_by_stable_id( "ENST00000217347" );
$tr->edits_enabled(1);

$tr->add_Attributes(
  Bio::EnsEMBL::Attribute->new(
    -code => '_rna_edit',
    -value => "1 0 CGTCGATGTTG",
  )
);
is(substr($tr->translateable_seq, 0, 6), 'ATGGCA', 'translation start in a seq_edit - seq before');
is(length($tr->translateable_seq),       804,      'translation start in a seq_edit - length before');

$tr->add_Attributes(
  Bio::EnsEMBL::Attribute->new(
    -code => '_transl_start',
    -value => "6",
  )
);
is(substr($tr->translateable_seq, 0, 6), 'ATGTTG', 'translation start in a seq_edit - seq after');
is(length($tr->translateable_seq),       874,      'translation start in a seq_edit - length after');

$tr->add_Attributes(
  Bio::EnsEMBL::Attribute->new(
    -code => '_rna_edit',
    -value => "869 868 CGTCGTGATTG",
  )
);
is(substr(reverse($tr->translateable_seq), 0, 11), reverse('CGTCGTGATTG'), 'translation end in a seq_edit - seq before');
is(length($tr->translateable_seq),                 885,                    'translation end in a seq_edit - length before');

$tr->add_Attributes(
  Bio::EnsEMBL::Attribute->new(
    -code => '_transl_end',
    -value => "884",
  )
);
is(substr(reverse($tr->translateable_seq), 0, 11), reverse('GAGTATCGTCG'), 'translation end in a seq_edit - seq after');
is(length($tr->translateable_seq),                 879,                    'translation end in a seq_edit - length after');

$multi->restore('core');

$tr = $ta->fetch_by_stable_id( "ENST00000217347" );
$tr->edits_enabled(1);

is(length($tr->translateable_seq), 804, 'explicit translation end with no seq_edit - length before');

$tr->add_Attributes(
  Bio::EnsEMBL::Attribute->new(
    -code => '_transl_end',
    -value => "873",
  )
);
is(length($tr->translateable_seq), 804, 'explicit translation end with no seq_edit - length after');

$tr->add_Attributes(
  Bio::EnsEMBL::Attribute->new(
    -code => '_rna_edit',
    -value => "869 868 CGTCGTGATTG",
  )
);
is(length($tr->translateable_seq), 809, 'explicit translation end with seq_edit - length after');

$multi->restore('core');

#
# tests for multiple versions of transcripts in a database
#

$tr = $ta->fetch_by_stable_id('ENST00000355555');
is( $tr->dbID, 21740, 'Fetched transcript by stable id' );

@transcripts = @{ $ta->fetch_all_versions_by_stable_id('ENST00000355555') };
is( scalar(@transcripts), 1, 'Fetched all transcripts by stable_id' );

$tr = $ta->fetch_by_translation_stable_id('ENSP00000355555');
is( $tr->dbID, 21740, 'Fetched transcript by translation stable id' );

$tr->stable_id_version('ENSP00000171455.4');
is($tr->stable_id, 'ENSP00000171455', 'Stable id set with stable_id_version');
is($tr->version, 4, 'Version set with stable_id_version');
is($tr->stable_id_version, 'ENSP00000171455.4', 'Stable id and version from stable_id_version');

$tr->stable_id_version('ENSP00000171456');
is($tr->stable_id, 'ENSP00000171456', 'Stable id set with stable_id_version');
is($tr->version, undef, 'Version undef from stable_id_version');
is($tr->stable_id_version, 'ENSP00000171456', 'Stable id and no version from stable_id_version');

$tr = $ta->fetch_by_translation_stable_id('ENSP00000355555.1');
is( $tr->dbID, 21740, 'Fetched transcript by translation stable id with version' );

$tr = $ta->fetch_by_translation_stable_id('ENSP00000355555.1a');
ok( ! defined($tr), 'Fetched transcript by translation stable id with bad version' );

$tr = $ta->fetch_by_translation_stable_id_version('ENSP00000355555', 1);
is( $tr->dbID, 21740, 'Fetched transcript by translation stable id (version) with version' );

$tr = $ta->fetch_by_translation_stable_id_version('ENSP00000355555', '1a');
ok( ! defined($tr), 'Fetched transcript by translation stable id (version) with bad version' );

@transcripts = @{ $ta->fetch_all_by_exon_stable_id('ENSE00001109603') };
is( scalar(@transcripts), 1);
is( $transcripts[0]->dbID, 21740, 'Fetched transcript by exon stable id' );

$g = $db->get_GeneAdaptor->fetch_by_stable_id('ENSG00000355555');
@transcripts = @{ $ta->fetch_all_by_Gene($g) };
is( scalar(@transcripts), 1);
is( $transcripts[0]->dbID, 21740, 'Fetched transcript by gene' );

my $sl = $sa->fetch_by_region('chromosome', 'MT_NC_001807');
@transcripts = @{ $sl->get_all_Transcripts };
is( scalar(@transcripts), 1, 'Fetched all transcripts for region' );

@transcripts = @{ $ta->fetch_all_by_external_name('MAE1_HUMAN') };
is( scalar(@transcripts), 1);
is( $transcripts[0]->dbID, 21738, 'Fetched transcript by external name' );

$tr = $ta->fetch_by_display_label('MAPRE1');
is( $tr->dbID, 21738, 'Fetched transcript by display label' );

# store/update

$tr = $ta->fetch_by_stable_id('ENST00000355555');
$g = $db->get_GeneAdaptor->fetch_by_transcript_id($tr->dbID);
$tr->get_all_Exons;

my $tl = $tr->translation;

$multi->hide( "core", "gene", "transcript", "translation", "meta_coord" );

$tr->version(3);
$tr->dbID(undef);
$tr->adaptor(undef);
$ta->store($tr, $g->dbID);

$tr->version(4);
$tr->is_current(0);
$tr->dbID(undef);
$tr->adaptor(undef);
$ta->store($tr, $g->dbID);

$tr->version(undef);
$tr->is_current(0);
$tr->dbID(undef);
$tr->adaptor(undef);
$tl->version(undef);
$tr->translation($tl);
$ta->store($tr, $g->dbID);

$tr = $ta->fetch_by_stable_id('ENST00000355555');
is($tr->is_current, 1, 'Transcript is current');   # 148

@transcripts = @{ $ta->fetch_all_versions_by_stable_id('ENST00000355555') };
foreach my $t (@transcripts) {
  if (defined $t->version && $t->version == 4) {
    is($t->is_current, 0, 'Transcript is not current');  # 149
  }
}

$tr->is_current(0);
$ta->update($tr);
my $t1 = $ta->fetch_by_stable_id('ENST00000355555');
ok(!$t1);   # 150

$tr->is_current(1);
$ta->update($tr);
$tr = $ta->fetch_by_stable_id('ENST00000355555');
is($tr->is_current, 1, 'Transcript is now current');   # 151

my $null_versions = 0;
foreach my $t (@transcripts) {
  if (! defined $t->version) {
    ok(! defined $t->translation->version);
    $null_versions++;
  }
}
is ( $null_versions, 1, "Null/undef version stored and retrieved");

$multi->restore;

# UTR Tests
{
  my $utr_testing = sub {
    my ($tid, $five_start, $five_end, $three_start, $three_end) = @_;
    my $t;
    if(ref($tid)) {
      $t = $tid;
    }
    else {
      $t = $db->get_TranscriptAdaptor()->fetch_by_dbID($tid);
    }
    my $five_utr = $t->five_prime_utr_Feature();
    my $three_utr = $t->three_prime_utr_Feature();
    is($five_utr->start(), $five_start, 'Checking five prime UTR start');
    is($five_utr->end(), $five_end, 'Checking five prime UTR end');
    is($three_utr->start(), $three_start, 'Checking three prime UTR start');
    is($three_utr->end(), $three_end, 'Checking three prime UTR end');
  };
  my $three_prime_seq_test = sub {
    my ($tid) = @_;
    my $t = $db->get_TranscriptAdaptor()->fetch_by_dbID($tid);
    my $true_seq = $t->three_prime_utr()->seq();
    my $seq = $t->three_prime_utr_Feature()->feature_Slice()->seq();
    is($seq, $true_seq, 'Asserting 3 prime seq as expected; we have a transcript ending in the last Exon');
  };
  note 'Testing +ve strand 5 prime Exon UTR readthrough';
  $utr_testing->(21729, 30653524, 30685637, 30707177, 30709209);
  
  note 'Testing -ve strand 5 prime Exon UTR readthrough';
  $utr_testing->(21726, 30578039, 30583588, 30568364, 30572314);
  $three_prime_seq_test->(21726);

  my $tid = 21729;
  my $t = $db->get_TranscriptAdaptor()->fetch_by_dbID($tid);
  my $five_utrs = $t->get_all_five_prime_UTRs();
  is(scalar(@$five_utrs), 2, "There are 2 five prime utr features");
  is($five_utrs->[0]->start(), $t->seq_region_start(), "Correct five prime UTR start");
  is($five_utrs->[1]->end(), 30685637, "Correct five prime UTR end");

  my $three_utrs = $t->get_all_three_prime_UTRs();
  is(scalar(@$three_utrs), 1, "There is one three prime utr feature");
  is($three_utrs->[0]->start(), 30707177, "Correct three prime UTR start");
  is($three_utrs->[0]->end(), $t->seq_region_end(), "Correct three prime UTR end");

  $tid = 21726;
  $t = $db->get_TranscriptAdaptor()->fetch_by_dbID($tid);
  $five_utrs = $t->get_all_five_prime_UTRs();
  is(scalar(@$five_utrs), 2, "There are 2 five prime utr features on reverse strand");
  is($five_utrs->[1]->start(), 30578039, "Correct five prime UTR start on reverse strand");
  is($five_utrs->[0]->end(), $t->seq_region_end(), "Correct five prime UTR end on reverse strand");

  $three_utrs = $t->get_all_three_prime_UTRs();
  is(scalar(@$three_utrs), 1, "There is one three prime utr feature on reverse strand");
  is($three_utrs->[0]->start(), $t->seq_region_start(), "Correct three prime UTR start on reverse strand");
  is($three_utrs->[0]->end(), 30572314, "Correct three prime UTR end on reverse strand");
  
  # we have to build some of our own as it's easier to see the coordinates & do the maths
  my $transcript_builder = sub {
    my ($local_slice, $exon_array, $seq_start, $seq_end, $exon_start_idx, $exon_end_idx) = @_;
    my @exons = map {
      Bio::EnsEMBL::Exon->new(-START => $_->[0], -END => $_->[1], -SLICE => $local_slice, -STRAND => $_->[2]);
    } @{$exon_array};
    my $t = Bio::EnsEMBL::Transcript->new(-SLICE => $local_slice, -EXONS => \@exons);
    $t->translation(Bio::EnsEMBL::Translation->new(
      -SEQ_START => $seq_start,
      -SEQ_END => $seq_end,
      -START_EXON => $exons[$exon_start_idx],
      -END_EXON => $exons[$exon_end_idx],
    ));
    return $t;
  };
  
  #Local slice of the whole of Chr20
  my $local_slice = $sa->fetch_by_region('chromosome', "20");
  
  note 'Testing +ve strand 3 prime Exon UTR readthrough';
  my $read_through_pos_three_prime = $transcript_builder->($local_slice, [[1,9,1],[30,39,1],[90,99,1]], 2, 3, 0, 1);
  $utr_testing->($read_through_pos_three_prime, 1, 1, 33, 99);

  note 'Testing -ve strand prime Exon UTR readthrough';
  my $read_through_neg_three_prime = $transcript_builder->($local_slice, [[90,99,-1], [30,39,-1], [1,9,-1]], 2, 3, 0, 1);
  $utr_testing->($read_through_neg_three_prime, 99, 99, 1, 36);

  note 'Testing no UTRs; we expect no features or sequence to be returned';
  my $no_pos_utrs = $transcript_builder->($local_slice, [[1,9,1]], 1, 9, 0, 0);
  ok(! defined $no_pos_utrs->five_prime_utr_Feature(), 'No 5 prime UTR means no feature');
  ok(! $no_pos_utrs->five_prime_utr(), 'No 5 prime UTR means no seq');
  ok(! defined $no_pos_utrs->three_prime_utr_Feature(), 'No 3 prime UTR means no feature');
  ok(! defined $no_pos_utrs->three_prime_utr(), 'No 3 prime UTR means no seq');
  
  note 'Testing no UTRs; we expect no features or sequence to be returned';
  my $no_pos_neg_utrs = $transcript_builder->($local_slice, [[1,9,-1]], 1, 9, 0, 0);
  ok(! defined $no_pos_neg_utrs->five_prime_utr_Feature(), 'No 5 prime UTR means no feature');
  ok(! $no_pos_neg_utrs->five_prime_utr(), 'No 5 prime UTR means no seq');
  ok(! defined $no_pos_neg_utrs->three_prime_utr_Feature(), 'No 3 prime UTR means no feature');
  ok(! defined $no_pos_neg_utrs->three_prime_utr(), 'No 3 prime UTR means no seq');
  
  
}

# CDS tests

{

  my $tid = 21726;
  my $t = $db->get_TranscriptAdaptor()->fetch_by_dbID($tid);
  my $cds = $t->get_all_CDS();
  is(scalar(@$cds), scalar(@{$t->get_all_translateable_Exons()}), "There are 2 coding structures");
  is($cds->[1]->start, $t->coding_region_start, "Correct coding start");
  is($cds->[0]->end, $t->coding_region_end, "Correct coding end");

}

SKIP: {
  skip 'No registry support for SQLite yet', 1 if $db->dbc->driver() eq 'SQLite';

  #test the get_species_and_object_type method from the Registry
  my $registry = 'Bio::EnsEMBL::Registry';
  my ( $species, $object_type, $db_type ) = $registry->get_species_and_object_type('ENST00000355555');
  ok( $species eq 'homo_sapiens' && $object_type eq 'Transcript');
}


## Relative vs absolute coordinates test
print "Comparing relative and absolute coordinates\n";
## Retrieve transcript by id
my $tid = 21726;
my $absolute_transcript = $db->get_TranscriptAdaptor()->fetch_by_dbID($tid);
my @absolute_coords = sort { $a->end() <=> $b->end() } ( $absolute_transcript->genomic2pep($absolute_transcript->seq_region_start, $absolute_transcript->seq_region_end, $absolute_transcript->strand) );
## Retrieve same transcript via feature slice
my $relative_slice = $absolute_transcript->feature_Slice();
my $relative_transcripts = $relative_slice->get_all_Transcripts();
my $relative_transcript;
foreach my $transcript (@$relative_transcripts) {
  if ($transcript->stable_id eq $absolute_transcript->stable_id) {
    $relative_transcript = $transcript;
    last;
  }
}

my @relative_coords = sort { $a->end() <=> $b->end() } ( $relative_transcript->genomic2pep($relative_transcript->seq_region_start, $relative_transcript->seq_region_end, $relative_transcript->strand) );
is(scalar(@absolute_coords), scalar(@relative_coords), "Same number of results");

## Compare coordinates of mappings
for (my $i = 0; $i < scalar(@absolute_coords); $i++) {
  if ($absolute_coords[$i]->isa('Bio::EnsEMBL::Mapper::Gap')) {
    is(ref($relative_coords[$i]), 'Bio::EnsEMBL::Mapper::Gap', "Both are gaps");
  } else {
    is($absolute_coords[$i]->start, $relative_coords[$i]->start, "Starts match");
    is($absolute_coords[$i]->end, $relative_coords[$i]->end, "Ends match");
    is($absolute_coords[$i]->strand, $relative_coords[$i]->strand, "Strands match");
  }
}

done_testing();

#
# end main
#

sub test_trans_mapper_edits {
  $tr->edits_enabled(1);

  # We want to fetch the first exon, and edit it's coordinates
  # to test the mapper and edits
  my $exons = $tr->get_all_Exons();
  ok(@{$exons}, 'We found at least one exon in the transcript');
  my $first_exon = shift @{$exons};

  # Start and end for testing genomic to cda coordinates should be relative
  # to the first exon, whichever direction transcription is occurring
  my $start = ($tr->strand() == 1) ? $first_exon->seq_region_start() : $first_exon->seq_region_end() - 11;
  my $end = ($tr->strand() == 1) ? $first_exon->seq_region_start() + 11 : $first_exon->seq_region_end();

  @coords = $tr->genomic2cdna($start, $end, $tr->strand());

  ok(@coords == 1 && $coords[0]->isa('Bio::EnsEMBL::Mapper::Coordinate'));
  is($coords[0]->start(), 1, 'Start is correct');
  is($coords[0]->end(), 12, 'End is correct');

  # deletion of 3 bp
  my $se = Bio::EnsEMBL::SeqEdit->new
    (-CODE  => '_rna_edit',
     -START => 2,
     -END   => 4,
     -ALT_SEQ => '');

  $tr->add_Attributes($se->get_Attribute());

  @coords = $tr->genomic2cdna($start, $end, $tr->strand());

  note("Expect coord, gap, coord");
  print_coords(\@coords);

  ok(@coords == 3);
  ok($coords[0]->start == 1);
  ok($coords[0]->end   == 1);

  ok($coords[1]->isa('Bio::EnsEMBL::Mapper::Gap'));
  ok($coords[1]->length() == 3);

  ok($coords[2]->start == 2);
  ok($coords[2]->end   == 9);


  # replacement, should have no effect
  $se->start(6);
  $se->end(8);
  $se->alt_seq('ACT');
  $tr->add_Attributes($se->get_Attribute());

  @coords = $tr->genomic2cdna($start, $end, $tr->strand());

  note("Expect coord, gap, coord");
  print_coords(\@coords);

  ok(@coords == 3);
  ok($coords[0]->start == 1);
  ok($coords[0]->end   == 1);

  ok($coords[1]->isa('Bio::EnsEMBL::Mapper::Gap'));
  ok($coords[1]->length() == 3);

  ok($coords[2]->start == 2);
  ok($coords[2]->end   == 9);

  # insertion in middle

  $se->start(10);
  $se->end(9);
  $se->alt_seq('GGGG');
  $tr->add_Attributes($se->get_Attribute());

  @coords = $tr->genomic2cdna($start, $end, $tr->strand());

  note("Expect coord, gap, coord, coord");

  ok(@coords == 4);

  ok($coords[0]->start == 1);
  ok($coords[0]->end   == 1);

  ok($coords[1]->isa('Bio::EnsEMBL::Mapper::Gap'));
  ok($coords[1]->length() == 3);

  ok($coords[2]->start == 2);
  ok($coords[2]->end   == 6);

  ok($coords[3]->start == 11);
  ok($coords[3]->end   == 13);

  print_coords(\@coords);

  # insert at very start of cdna

  $se->start(1);
  $se->end(0);
  $se->alt_seq('A');
  $tr->add_Attributes($se->get_Attribute());

  @coords = $tr->genomic2cdna($start, $end, $tr->strand());

  note("Expect coords, gap, coord, coord");

  ok(@coords == 4);

  ok($coords[0]->start == 2);
  ok($coords[0]->end   == 2);

  ok($coords[1]->isa('Bio::EnsEMBL::Mapper::Gap'));
  ok($coords[1]->length() == 3);

  ok($coords[2]->start == 3);
  ok($coords[2]->end   == 7);

  ok($coords[3]->start == 12);
  ok($coords[3]->end   == 14);

  print_coords(\@coords);
}


#sleep 10;
# test fetching by supporting evidence
#$transcripts = $transcript_adaptor->fetch_all_by_exon_supporting_evidence(hit,feat);
#$transcripts = $transcript_adaptor->fetch_all_by_exon_supporting_evidence(hit,feat,anal);

sub print_coords {
  my $coords = shift;

  foreach my $c (@$coords) {
    if($c->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
      note("COORD ",$c->start .'-'.$c->end);
    } else {
      note("GAP (". $c->length().")");
    }
  }
}


