use lib 't';
use strict;
use warnings;

BEGIN { $| = 1;
	use Test;
	plan tests => 58;
}

use MultiTestDB;
use TestUtils qw ( debug test_getter_setter count_rows);
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::DnaDnaAlignFeature;

# switch on the debug prints
our $verbose = 0;

debug( "Startup test" );
ok(1);

my $multi = MultiTestDB->new();

my $db = $multi->get_DBAdaptor( "core" );

debug( "Test database instatiated" );
ok( $db );

my $gene;
my $ga = $db->get_GeneAdaptor();

debug ("Gene->list_dbIDs");
my $ids = $ga->list_dbIDs();
ok (@{$ids});

debug ("Gene->list_stable_ids");
my $stable_ids = $ga->list_stable_ids();
ok (@{$stable_ids});

$gene = $ga->fetch_by_stable_id( "ENSG00000171456" );

debug( "Gene->fetch_by_stable_id()" );
ok( $gene );



debug( "Gene dbID: ". $gene->dbID());
ok( $gene->dbID() == 18267 );

debug( "Gene start: ".$gene->start );
ok( $gene->start() == 30735607 );

debug( "Gene end: ".$gene->end );
ok( $gene->end() == 30815178 );

debug( "Gene external name: " . $gene->external_name );
ok( $gene->external_name eq "Q9H466");

debug( "Gene external dbname: " . $gene->external_db );
ok( $gene->external_db eq "SPTREMBL");

debug( "Gene display xref id: " . $gene->display_xref->dbID );
ok( $gene->display_xref->dbID() == 128324);


# test the getters and setters
ok( test_getter_setter( $gene, "external_name", "banana" ));   
ok( test_getter_setter( $gene, "external_db", "dummy" ));   
ok( test_getter_setter( $gene, "display_xref", 42 ));   


my $links = $gene->get_all_DBLinks();
debug( "Links: ".scalar( @$links ));

ok( scalar @$links == 6 );

# now create a new gene ...


my $sa = $db->get_SliceAdaptor();

my $slice = $sa->fetch_by_region( "chromosome", "20", 30_249_935, 31_254_640 );

debug( "Slice from SliceAdaptor" );
ok($slice);

my $analysis = $db->get_AnalysisAdaptor->fetch_by_logic_name("ensembl");
my $f_analysis = $db->get_AnalysisAdaptor->fetch_by_logic_name("Vertrna");
debug( "Analysis from AnalysisAdaptor" );
ok($analysis);


$gene = Bio::EnsEMBL::Gene->new();

my $transcript1 = Bio::EnsEMBL::Transcript->new();
my $transcript2 = Bio::EnsEMBL::Transcript->new();

my $ex1 = Bio::EnsEMBL::Exon->new(); 
my $ex2 = Bio::EnsEMBL::Exon->new();
my $ex3 = Bio::EnsEMBL::Exon->new();

my $translation1 = Bio::EnsEMBL::Translation->new();
my $translation2 = Bio::EnsEMBL::Translation->new();	

ok($gene);


$ex1->start( 13586 );
$ex1->end( 13735 );
$ex1->phase(0);
$ex1->end_phase( 0 );
$ex1->slice( $slice );
$ex1->strand(1);
$ex1->analysis($analysis);

my @feats;
my $fp = new Bio::EnsEMBL::FeaturePair;

$fp->start(13586);
$fp->end  (13705);
$fp->strand(1);
$fp->score(10);
$fp->slice($slice);
$fp->hstart(100);
$fp->hend    (219);
$fp->hstrand (1);
$fp->hseqname('dummy-hid');

push(@feats,$fp);


$fp = new Bio::EnsEMBL::FeaturePair;
$fp->start(13707);
$fp->end  (13735);
$fp->strand(1);
$fp->score(10);
$fp->slice($slice);

$fp->hstart  (220);
$fp->hend    (248);
$fp->hstrand (1);
$fp->hseqname('dummy-hid');
push(@feats,$fp);

#
#
# 2 Test DnaDnaAlignFeature::new(-features)
#
my $dnaf = Bio::EnsEMBL::DnaDnaAlignFeature->new( -features => \@feats );
$dnaf->analysis( $f_analysis );

$ex1->add_supporting_features( $dnaf );



$ex2->start(201372);
$ex2->end(201571);
$ex2->phase(0);
$ex2->end_phase( -1 );
$ex2->slice( $slice );
$ex2->strand(1);
$ex2->analysis($analysis);

@feats = ();
$fp = new Bio::EnsEMBL::FeaturePair;

$fp->start(201372);
$fp->end  (201471);
$fp->strand(1);
$fp->score(10);
$fp->slice( $slice );
$fp->hstart(100);
$fp->hend    (199);
$fp->hstrand (1);
$fp->hseqname('dummy-hid');

push(@feats,$fp);


$fp = new Bio::EnsEMBL::FeaturePair;
$fp->start(201472);
$fp->end  (201571);
$fp->strand(1);
$fp->score(10);
$fp->slice( $slice );


$fp->hstart  (201);
$fp->hend    (300);
$fp->hstrand (1);
$fp->hseqname('dummy-hid');
push(@feats,$fp);

#
#
# 2 Test DnaDnaAlignFeature::new(-features)
#
$dnaf = Bio::EnsEMBL::DnaDnaAlignFeature->new( -features => \@feats );
$dnaf->analysis( $f_analysis );

$ex2->add_supporting_features( $dnaf );

$ex3->start(210600);
$ex3->end(210800);
$ex3->phase(-1);
$ex3->end_phase( -1 );
$ex3->slice( $slice );
$ex3->strand(1);
$ex3->analysis($analysis);

$transcript1->add_Exon($ex1);
$transcript1->add_Exon($ex2);
$translation1->start_Exon($ex1);
$translation1->end_Exon($ex2);
$translation1->start(1);
$translation1->end(150);
$transcript1->translation($translation1);


$transcript2->add_Exon($ex1);
$transcript2->add_Exon($ex2);
$transcript2->add_Exon($ex3);
$translation2->start_Exon($ex1);
$translation2->end_Exon($ex2);
$translation2->start(1);
$translation2->end(180);
$transcript2->translation($translation2);

debug( "Transcripts created" );
ok($transcript1);


$gene->add_Transcript($transcript1);
$gene->add_Transcript($transcript2);

$gene->analysis($analysis);

debug( "Getting all the Transcripts/Exons from new Gene" );

my $count = 0;
my $translates  = 1;

foreach my $tr( @{$gene->get_all_Transcripts()} ) {
  if( $tr->translate()->seq() =~ /\*./ ) {
    $translates = 0;
    debug( "Translate failed." );
  }
  debug( "Translation: ".$tr->translate()->seq() );
  foreach my $exon ( @{$tr->get_all_Exons()} ) {
    debug( "  Exon start: ". $exon->start());
    debug( "  Exon end:   ". $exon->end() );
    debug( "  Exon strand ".$exon->strand() );
    $count++;
  }	
}

ok($count == 5);
ok( $translates );

ok( scalar(@{$gene->get_all_Exons()} ) == 3);

$gene = $gene->transform( "chromosome" );

my $desc = 'test description for a gene';
my $stable_id = 'ENSG00000171456';
$gene->description($desc);
$gene->stable_id($stable_id);

$multi->hide( "core", "meta_coord", "gene", "transcript", "exon", "exon_transcript", "gene_description", "translation", "gene_stable_id", "transcript_stable_id", "exon_stable_id", "translation_stable_id", "supporting_feature", "dna_align_feature" );

my $gene_ad = $db->get_GeneAdaptor();
debug( "Storing the gene" );
$gene_ad->store($gene);

ok(1);

my $genes = $slice->get_all_Genes();


ok(scalar( @$genes) == 1 );

my $gene_out = $genes->[0];

#make sure the stable_id was stored
ok($gene_out->stable_id eq $stable_id);

#make sure the description was stored
ok($gene_out->description eq $desc);

ok(scalar(@{$gene_out->get_all_Exons()}) == 3);


foreach my $tr( @{$gene_out->get_all_Transcripts()} ) {
  debug( "NewTranscript: ".$tr->dbID() );
  foreach my $exon ( @{$tr->get_all_Exons()} ) {
    debug( "  NewExon: ".$exon->start(). " ".$exon->end()." ".$exon->strand());
  }	
}

my $exons = $gene_out->get_all_Transcripts()->[0]->get_all_Exons();

ok( $exons->[0]->start == 13586 );


ok( $exons->[1]->strand == 1 );
ok( $exons->[1]->phase == 0 );



my $pep;
my $translate = 0;
foreach my $trans( @{$gene_out->get_all_Transcripts()} ){

  my $pep = $trans->translate();
  debug( "Peptide: ".$pep->seq() );

  if($pep->seq !~ /\*./){
    $translate = 1;
  } else {
    $translate = 0;
  }
}

ok($translate == 1);

my $t = $gene_out->get_all_Transcripts()->[1];

my $e = $t->get_all_Exons()->[0];
my $se = $e->get_all_supporting_features();

debug( "Got ".scalar( @$se )." supporting features." );
ok( scalar( @$se ) == 1 );

my $se_start = $se->[0]->start(); 


my $se_end = $se->[0]->end();


debug( "Supporting start $se_start, end $se_end" );
debug( "Exon start ".$e->start()." end ".$e->end() );

ok( $se_start == $e->start() );
ok( $se_end == $e->end() );


my $pep1 = $t->translate()->seq();

$e->phase(1);
my $pep2 = $t->translate()->seq();

debug( "Pep phase 0: $pep1" );
debug( "Pep phase 1: $pep2" );

ok( $pep1 ne $pep2 );
debug( "checking external references" );

$multi->restore();

$slice = $db->get_SliceAdaptor()->fetch_by_region
  ( "chromosome", "20", 30_252_000, 31_252_001 );

my $known = 0;
my $unknown = 0;

$genes = $slice->get_all_Genes();
for my $gene ( @$genes ) {
  if( $gene->is_known() ) {
    $known++;
  } else {
    $unknown++;
  }
}

debug( "known: $known Unknown: $unknown\n" );

ok( $known==17 );

#if( my $lite = $multi->get_DBAdaptor( 'lite' ) ) {
#  debug( "Lite database available" );
#  my $lga = $lite->get_GeneAdaptor();
#  $gene = $ga->fetch_by_stable_id( "ENSG00000171456" );

#  $lga->store( $gene );
#  debug( "Store done" );
#}


#save contents of gene table
$multi->save('core', 'gene');

# tests for update method
# go get a fresh gene again
$gene = $ga->fetch_by_stable_id( "ENSG00000171456" ); 

# the first update should no effect
$ga->update($gene);

my $newgene = $ga->fetch_by_stable_id( "ENSG00000171456" ); 
ok ( $newgene->display_xref->dbID() == 128324 );
ok ( $newgene->type eq 'ensembl' );

# now change the original gene and update it
my $dbEntryAdaptor=  $db->get_DBEntryAdaptor();

$gene->display_xref( $dbEntryAdaptor->fetch_by_dbID( 614 ));
$gene->type('dummy');
$ga->update($gene);

$newgene = $ga->fetch_by_stable_id( "ENSG00000171456" );
ok ( $newgene->display_xref->dbID() == 614 );
ok ( $newgene->type eq 'dummy' );

$multi->restore('core', 'gene');


#
# test GeneAdaptor::fetch_all_by_domain
#
my @genes = @{$ga->fetch_all_by_domain('IPR000010')};

debug("Fetch by domain 'IPR000010'");

ok(@genes == 2);
debug("Got " . scalar(@genes) . " genes");
ok(($genes[0]->stable_id() eq 'ENSG00000131044') ||
   ($genes[1]->stable_id() eq 'ENSG00000131044'));
ok(($genes[0]->stable_id() eq 'ENSG00000174873') ||
   ($genes[1]->stable_id() eq 'ENSG00000174873'));


#
# test GeneAdaptor::fetch_all_by_external_name
#

#Q15691
($gene) = @{$ga->fetch_all_by_external_name('MAE1_HUMAN')};
debug($gene->stable_id);
ok($gene->stable_id() eq 'ENSG00000101367');

#
# test GeneAdaptor::get_Interpro_by_geneid
#
debug("Test get_Interpro_by_geneid");
my @interpro = @{$ga->get_Interpro_by_geneid('ENSG00000174873')};
ok(@interpro == 1);
debug($interpro[0]);


ok($gene->display_id eq $gene->stable_id);

#
# test Gene: get_all_alt_alleles
#

$gene = $ga->fetch_by_dbID( 18256 );
my $alt_genes = $gene->get_all_alt_alleles();

ok( scalar( @$alt_genes ) == 3 );

# expect the following alleles
my %gene_ids = ( 18257 => 1, 18258 => 1, 18259 => 1);
my $ok = 1;
for my $gene ( @$alt_genes ) {
  $ok = $ok && $gene_ids{$gene->dbID()};
}
ok( $ok );

#
# test storing a new allele group
#

$multi->hide( 'core', 'alt_allele' );

my @alt_genes;
push( @alt_genes, $ga->fetch_by_dbID(18270) );
push( @alt_genes, $ga->fetch_by_dbID(18271) );
push( @alt_genes, $ga->fetch_by_dbID(18272) );
$ga->store_alt_alleles( \@alt_genes ); 

$gene = $ga->fetch_by_dbID( 18270 );
$alt_genes = $gene->get_all_alt_alleles();
%gene_ids = ( 18271=>1, 18272=>1 );

$ok = 1;
for my $gene ( @$alt_genes ) {
  $ok = $ok && $gene_ids{$gene->dbID()};
}
ok( $ok );

#
# Gene remove test
#

$multi->save( "core", "gene", "gene_stable_id", "gene_description",
	      "transcript", "transcript_stable_id",
	      "translation", "translation_stable_id", "protein_feature",
	      "exon", "exon_stable_id", "exon_transcript", "supporting_feature",
	      "object_xref", "go_xref", "identity_xref" );

$gene = $ga->fetch_by_stable_id( "ENSG00000171456" );

my $gene_count = count_rows( $db, "gene" );
my $exon_count = count_rows( $db, "exon" );
my $trans_count = count_rows( $db, "transcript" );
my $tl_count = count_rows( $db, "translation" );
my $gdesc_count = count_rows($db, "gene_description");
my $gstable_count = count_rows($db, "gene_stable_id");

my $tminus = scalar( @{$gene->get_all_Transcripts() } );
my $eminus = scalar( @{$gene->get_all_Exons() } );

debug( "Genes before ".$gene_count );
debug( "Exons before ".$exon_count );
debug( "Transcripts before ".$trans_count );
debug( "Translations before ".$tl_count );
debug( "Gene has ".$tminus." transcripts" );
debug( "Gene has ".$eminus." exons" );

$ga->remove( $gene );

ok( count_rows( $db, "gene" ) == ( $gene_count - 1 ));
ok( count_rows( $db, "transcript" ) == ($trans_count-$tminus));
ok( count_rows( $db, "exon" ) == ( $exon_count - $eminus ));
ok( count_rows( $db, "gene_description" ) == ( $gdesc_count -1 ));
ok( count_rows( $db, "gene_stable_id" ) == ($gstable_count -1));

ok(!defined($gene->dbID()));
ok(!defined($gene->adaptor()));

$multi->restore('core');

#
# regression test - test the recalculation of coords
# in the Gene.  This was setting the end incorrectly
# before.
#
$gene = Bio::EnsEMBL::Gene->new();

$gene->slice($slice);

my $first_ex = Bio::EnsEMBL::Exon->new
  (-START => 10,
   -END   => 100,
   -STRAND => 1,
   -SLICE => $slice);

my $second_ex = Bio::EnsEMBL::Exon->new
  (-START   => 200,
   -END     => 400,
   -STRAND  => 1,
   -SLICE   => $slice);

$transcript1 = Bio::EnsEMBL::Transcript->new
  (-EXONS => [$first_ex, $second_ex]);

$transcript2 = Bio::EnsEMBL::Transcript->new
  (-EXONS => [$first_ex]);


$gene->add_Transcript($transcript1);
$gene->add_Transcript($transcript2);

$gene->recalculate_coordinates();

ok($gene->start() == 10);
ok($gene->end()  == 400);
