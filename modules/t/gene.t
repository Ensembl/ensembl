use lib 't';
use strict;
use warnings;

BEGIN { $| = 1;  
	use Test;
	plan tests => 39;
}

use MultiTestDB;
use TestUtils qw ( debug test_getter_setter );

use Bio::EnsEMBL::Gene;

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

my $e_slice = Bio::EnsEMBL::Slice->new
  ( -empty => 1,
    -adaptor => $db->get_SliceAdaptor() 
  );

$gene->transform( $e_slice );

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

my $slice = $sa->fetch_by_chr_start_end( "20", 30_249_935, 31_254_640 );

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
$ex1->contig( $slice );
$ex1->strand(1);
$ex1->analysis($analysis);

my @feats;
my $fp = new Bio::EnsEMBL::FeaturePair;

$fp->start(13586);
$fp->end  (13705);
$fp->strand(1);
$fp->score(10);
$fp->contig($slice);
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
$fp->contig($slice);
$fp->seqname(1);

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
$ex2->contig( $slice );
$ex2->strand(1);
$ex2->analysis($analysis);

@feats = ();
$fp = new Bio::EnsEMBL::FeaturePair;

$fp->start(201372);
$fp->end  (201471);
$fp->strand(1);
$fp->score(10);
$fp->contig($slice);
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
$fp->contig($slice);
$fp->seqname(1);

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
$ex3->contig( $slice );
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
    debug( "Translate failed.".$tr->translate()->seq() );
  }

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

$gene->transform();
$multi->hide( "core", "gene", "transcript", "exon", "exon_transcript", "gene_description", "translation", "gene_stable_id", "transcript_stable_id", "exon_stable_id", "translation_stable_id", "supporting_feature", "dna_align_feature" );

my $gene_ad = $db->get_GeneAdaptor();
my $desc = 'test description';
$gene->description($desc);
debug( "Storing the gene" );

$gene_ad->store($gene);

ok(1);

my $genes = $slice->get_all_Genes();



ok(scalar( @$genes) == 1 );

my $gene_out = $genes->[0];

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

ok($gene_out->description eq $desc);


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
ok( scalar( @$se ) == 2 );

my $se_start = $se->[0]->start() < $se->[1]->start() ? 
  $se->[0]->start() : $se->[1]->start();

my $se_end = $se->[0]->end() > $se->[1]->end() ? 
  $se->[0]->end() : $se->[1]->end();


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

$slice = $db->get_SliceAdaptor()->fetch_by_chr_start_end("20", 30_252_000, 31_252_001 );

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
