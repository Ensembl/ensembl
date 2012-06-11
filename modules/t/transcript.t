use strict;
use warnings;

use Test::More;

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
ok($tr && $tr->dbID() == 21737 );


$tr = $ta->fetch_by_stable_id( "ENST00000217347" );

$tr = $tr->transform('contig');

ok( $tr );

note ( "External transcript name: " . $tr->external_name );
ok ( $tr->external_name eq "MAPRE1");

note ( "External transcript dbname: " . $tr->external_db );
ok ( $tr->external_db eq 'HUGO' );

note ( "Display_xref_id: " . $tr->display_xref->dbID() );
ok ( $tr->display_xref->dbID() == 97759 );
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

note( "start() == ".$tr->start() );
ok( $tr->start() == 79874 );

note( "end() == ".$tr->end() );
ok( $tr->end() == 110306 );

note( "spliced_seq->substr == \"".substr( $tr->spliced_seq(),0, 10 )."\"" );
ok( substr( $tr->spliced_seq(), 0, 10 ) eq "ACGAGACGAA" ); 

note( "translateable_seq->substr == \"".substr( $tr->translateable_seq(),0,10 )."\"" );
ok( substr( $tr->translateable_seq(),0,10 ) eq "ATGGCAGTGA" );

note( "coding_region_start() == ".$tr->coding_region_start() );
ok( $tr->coding_region_start() == 85834 );

note( "coding_region_end() == ".$tr->coding_region_end() );
ok( $tr->coding_region_end() == 108631 );

note( "pep2genomic: ".($tr->pep2genomic( 10,20 ))[0]->start());
my @pepcoords = $tr->pep2genomic( 10, 20 );
ok( $pepcoords[0]->start() == 85861 );

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

note("expecting $coord_num coords (in pep coords), got " . scalar(@coords));
ok(scalar(@coords) == $coord_num);
my ($last_coord) = reverse(@coords); 
note("expecting peptide length: $pep_len, got".$last_coord->end);
ok($pep_len == $last_coord->end);

note( "start Exon: ".$tr->start_Exon->stable_id() );
note( "end Exon: ".$tr->end_Exon->stable_id() );

note( "cdna_coding_start: ". $tr->cdna_coding_start );
ok($tr->cdna_coding_start == 65);
ok(test_getter_setter($tr, 'cdna_coding_start', 99));

note( "five_prime_utr: ".substr( $tr->five_prime_utr()->seq(), -5 , 5 ));
ok( substr( $tr->five_prime_utr()->seq(), -5, 5) eq "CGAAG" ); 

note( "cdna_coding_end: ". $tr->cdna_coding_end );
ok($tr->cdna_coding_end == 868);
ok(test_getter_setter($tr, 'cdna_coding_end', 102));

note( "three_prime_utr: ".substr( $tr->three_prime_utr()->seq(), -5, 5  ));
ok( substr( $tr->three_prime_utr()->seq(), -5, 5  ) eq "TTCAA");

note( "Transcript has: ". scalar( @{$tr->get_all_Exons()} ). " Exons" );
ok( scalar( @{$tr->get_all_Exons()} ) == 7 );

note( "Flushing Exons" );
$tr->flush_Exons();

ok( scalar( @{$tr->get_all_Exons()} ) == 0 );


# get a fresh tr to check the update method
$tr = $ta->fetch_by_stable_id( "ENST00000217347" );

$multi->save('core', 'transcript');

# the first update should have no effect
$ta->update($tr);

my $up_tr = $ta->fetch_by_stable_id( "ENST00000217347" );
ok ( $up_tr->display_xref->dbID() == 97759 );

my $dbentryAdaptor = $db->get_DBEntryAdaptor();

$tr->display_xref($dbentryAdaptor->fetch_by_dbID( 614 ));
$ta->update($tr);

$up_tr = $ta->fetch_by_stable_id( "ENST00000217347" );
ok ( $up_tr->display_xref->dbID() == 614 );

$multi->restore('core', 'transcript');



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
ok($tr && $tr->stable_id eq 'ENST00000246229');

#
# test fetch_by_translation_id
#

$tr = $ta->fetch_by_translation_id(21734);

ok($tr && $tr->stable_id eq 'ENST00000201961');

ok($tr->display_id() eq $tr->stable_id());

#
# test TranscriptAdaptor::fetch_all_by_biotype
#
note("Test fetch_all_by_biotype");
my @transcripts = @{$ta->fetch_all_by_biotype('protein_coding')};
ok(@transcripts == 25);
@transcripts = @{$ta->fetch_all_by_biotype(['protein_coding','pseudogene'])};
note "Got ".scalar(@transcripts)." transcripts\n";
ok(@transcripts == 25);


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

  my $orig_seq = $transcript->slice->subseq(
					    $transcript->start(),
					    $transcript->end(), 
					    $transcript->strand());

  my $idl=0;
  my $new_seq = $exons[0]->seq()->seq();
  foreach my $intron (@introns){
    $new_seq .= $intron->seq;
    $new_seq .= $exons[$idl+1]->seq->seq();
    $idl++;

  }

  ok($orig_seq eq $new_seq);

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
             "ontology_xref", "identity_xref");

$ta->remove($tr);

ok(!defined($tr->dbID()));
ok(!defined($tr->adaptor()));

ok( count_rows( $db, "transcript") == ($tr_count - 1));
ok( count_rows( $db, "translation") == ($tl_count - 1));
ok( count_rows( $db, "exon_transcript") == ($ex_tr_count - $ex_tr_minus));

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

my $seq2 = $tr->spliced_seq();

ok( $seq1 ne $seq2 );
ok( $seq2 =~ /^GATTACA/ );

my $cdna_cds_start2 = $tr->cdna_coding_start();
my $cdna_cds_end2   = $tr->cdna_coding_end();

ok($cdna_cds_start1 == $cdna_cds_start2 - 1);
ok($cdna_cds_end1   == $cdna_cds_end2   - 1);


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

ok($cdna_cds_start2 == $tr->cdna_coding_start());
ok($cdna_cds_end2   == $tr->cdna_coding_end() - 3);

# test that the edits can be disabled
$tr->edits_enabled(0);

ok($tr->cdna_coding_start() == $cdna_cds_start1);
ok($tr->cdna_coding_end()   == $cdna_cds_end1);

ok($tr->spliced_seq() eq $seq1);
ok($tr->translateable_seq() eq $tlseq1);

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

ok( $tr->translateable_seq() eq $tlseq2 );
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

ok($tr->translate->seq() eq 'MNFALILMINTLLALLLMIITFWLPQLNGYMEKSTPYECGFDPMSPARVPFSMKFFLVAITFLLFDLEIALLLPLPWALQTTNLPLMVMSSLLLIIILALSLAYEWLQKGLDWAE');



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
  ok( $tr->spliced_seq() eq $mapped_tr->spliced_seq() );
}

$multi->hide('core', 'transcript', 'transcript_attrib', 'translation',
             'exon_transcript', 'exon');


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

ok(count_rows($db, 'transcript_attrib') == 2);

$multi->restore('core');

#
# tests for multiple versions of transcripts in a database
#

$tr = $ta->fetch_by_stable_id('ENST00000355555');
note("fetch_by_stable_id");
ok( $tr->dbID == 21740 );

@transcripts = @{ $ta->fetch_all_versions_by_stable_id('ENST00000355555') };
note("fetch_all_versions_by_stable_id");
ok( scalar(@transcripts) == 1 );

$tr = $ta->fetch_by_translation_stable_id('ENSP00000355555');
note("fetch_by_translation_stable_id");
ok( $tr->dbID == 21740 );

@transcripts = @{ $ta->fetch_all_by_exon_stable_id('ENSE00001109603') };
note("fetch_all_by_exon_stable_id");
ok( scalar(@transcripts) == 1 && $transcripts[0]->dbID == 21740 );

$g = $db->get_GeneAdaptor->fetch_by_stable_id('ENSG00000355555');
@transcripts = @{ $ta->fetch_all_by_Gene($g) };
note("fetch_all_by_Gene");
ok( scalar(@transcripts) == 1 && $transcripts[0]->dbID == 21740 );

my $sl = $sa->fetch_by_region('chromosome', 'MT_NC_001807');
@transcripts = @{ $sl->get_all_Transcripts };
ok( scalar(@transcripts) == 1 );

@transcripts = @{ $ta->fetch_all_by_external_name('MAE1_HUMAN') };
note( "fetch_all_by_external_name" );
ok( scalar(@transcripts) == 1 && $transcripts[0]->dbID == 21738 );

$tr = $ta->fetch_by_display_label('MAPRE1');
note("fetch_by_display_label");
ok( $tr->dbID == 21738 );

# store/update

$tr = $ta->fetch_by_stable_id('ENST00000355555');
$g = $db->get_GeneAdaptor->fetch_by_transcript_id($tr->dbID);
$tr->get_all_Exons;

$multi->hide( "core", "gene", "transcript", "exon", 'xref', 'object_xref',
              "exon_transcript", "translation" );

$tr->version(3);
$tr->dbID(undef);
$tr->adaptor(undef);
$ta->store($tr, $g->dbID);

$tr->version(4);
$tr->is_current(0);
$tr->dbID(undef);
$tr->adaptor(undef);
$ta->store($tr, $g->dbID);
$tr = $ta->fetch_by_stable_id('ENST00000355555');
ok($tr->is_current == 1);   # 148

@transcripts = @{ $ta->fetch_all_versions_by_stable_id('ENST00000355555') };
foreach my $t (@transcripts) {
  next unless ($t->version == 4);
  ok($t->is_current == 0);  # 149
}

$tr->is_current(0);
$ta->update($tr);
my $t1 = $ta->fetch_by_stable_id('ENST00000355555');
ok(!$t1);   # 150

$tr->is_current(1);
$ta->update($tr);
$tr = $ta->fetch_by_stable_id('ENST00000355555');
ok($tr->is_current == 1);   # 151

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

{
  #test the get_species_and_object_type method from the Registry
  my $registry = 'Bio::EnsEMBL::Registry';
  my ( $species, $object_type, $db_type ) = $registry->get_species_and_object_type('ENST00000355555');
  ok( $species eq 'homo_sapiens' && $object_type eq 'Transcript');
}

done_testing();

#
# end main
#

sub test_trans_mapper_edits {
  $tr->edits_enabled(1);


  my $start = ($tr->strand() == 1) ? 1  : $tr->end() - 11;
  my $end   = ($tr->strand() == 1) ? 12 : $tr->end();

  @coords = $tr->genomic2cdna($start, $end, $tr->strand());

  ok(@coords == 1 && $coords[0]->isa('Bio::EnsEMBL::Mapper::Coordinate'));
  ok($coords[0]->start() == 1);
  ok($coords[0]->end()   == 12);

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


