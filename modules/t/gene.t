use lib 't';
use strict;
use warnings;

BEGIN { $| = 1;  
	use Test;
	plan tests => 22;
}

use MultiTestDB;
use TestUtils qw ( debug );

use Bio::EnsEMBL::Gene;

# switch on the debug prints
my $verbose = 0;

debug( "Startup test" );
ok(1);

my $multi = MultiTestDB->new();

my $db = $multi->get_DBAdaptor( "core" );

debug( "Test database instatiated" );
ok( $db );

my $gene;
my $ga = $db->get_GeneAdaptor();

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



my $links = $gene->get_all_DBLinks();
debug( "Links: ".scalar( @$links ));

ok( scalar @$links == 6 );

# now create a new gene ...

my $slice = $db->get_SliceAdaptor()->fetch_by_chr_start_end( "20",30264615, 30265615 );
debug( "Slice from SliceAdaptor" );
ok($slice);

my $analysis = $db->get_AnalysisAdaptor->fetch_by_logic_name("ensembl");
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


$ex1->start(5);
$ex1->end(10);
$ex1->phase(0);
$ex1->end_phase( 0 );
$ex1->contig( $slice );
$ex1->strand(1);
$ex1->analysis($analysis);

$ex2->start(15);
$ex2->end(23);
$ex2->phase(0);
$ex2->end_phase( 0 );
$ex2->contig( $slice );
$ex2->strand(1);
$ex2->analysis($analysis);

$ex3->start(28);
$ex3->end(33);
$ex3->phase(0);
$ex3->end_phase( 0 );
$ex3->contig( $slice );
$ex3->strand(1);
$ex3->analysis($analysis);

$transcript1->add_Exon($ex1);
$transcript1->add_Exon($ex2);
$translation1->start_Exon($ex1);
$translation1->end_Exon($ex2);
$translation1->start(1);
$translation1->end(9);
$transcript1->translation($translation1);


$transcript2->add_Exon($ex1);
$transcript2->add_Exon($ex2);
$transcript2->add_Exon($ex3);
$translation2->start_Exon($ex1);
$translation2->end_Exon($ex3);
$translation2->start(1);
$translation2->end(6);
$transcript2->translation($translation2);

debug( "Transcripts created" );
ok($transcript1);


$gene->add_Transcript($transcript1);
$gene->add_Transcript($transcript2);

$gene->analysis($analysis);

debug( "Getting all the Transcripts/Exons from new Gene" );

my $count = 0;

foreach my $tr( @{$gene->get_all_Transcripts()} ) {
  foreach my $exon ( @{$tr->get_all_Exons()} ) {
    debug( "  Exon start: ". $exon->start());
    debug( "  Exon end:   ". $exon->end() );
    debug( "  Exon strand ".$exon->strand() );
    $count++;
  }	
}

ok($count == 5);

ok( scalar(@{$gene->get_all_Exons()} ) == 3);

$gene->transform();
$multi->hide( "core", "gene", "transcript", "exon", "exon_transcript", "gene_description", "translation", "gene_stable_id", "transcript_stable_id", "exon_stable_id", "translation_stable_id" );

my $gene_ad = $db->get_GeneAdaptor();
debug( "Storing the gene" );
$gene_ad->store($gene);

ok(1);

my $new_slice = $db->get_SliceAdaptor()->fetch_by_chr_start_end( "20",30263615, 30265615 );

my $genes = $new_slice->get_all_Genes();



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

ok($exons->[0]->start==1005);
ok($exons->[1]->strand==1);
ok($exons->[1]->phase==0);



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

my $pep1 = $t->translate()->seq();

$e->phase(1);
my $pep2 = $t->translate()->seq();

debug( "Pep phase 0: $pep1" );
debug( "Pep phase 1: $pep2" );

ok( $pep1 ne $pep2 );
debug( "checking external references" );

$multi->restore();

$slice = $db->get_SliceAdaptor()->fetch_by_chr_start_end("20", 30_252_000, 31_252_001 );

my ( $known, $unknown );

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










