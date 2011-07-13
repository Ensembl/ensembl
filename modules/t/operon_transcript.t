use strict;
use warnings;

BEGIN {
	$| = 1;
	use Test;
	plan tests => 40;
}
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Operon;
use Bio::EnsEMBL::OperonTranscript;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::DBSQL::OperonAdaptor;

debug("Startup test");
ok(1);

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

my $dba = $multi->get_DBAdaptor("core");

debug("Test database instatiated");
ok($dba);

# get a slice
my $slice = $dba->get_SliceAdaptor()->fetch_by_seq_region_id(469283);

# create operon
my $start         = 3403162;
my $end           = 3405288;
my $strand        = 1;
my $display_label = "accBC";
my $operon = Bio::EnsEMBL::Operon->new( -START         => $start,
										-END           => $end,
										-STRAND        => $strand,
										-SLICE         => $slice,
										-DISPLAY_LABEL => $display_label, -STABLE_ID=>"op1" );

my $gene_name    = "accB";
my $gene_start   = 31225346;
my $gene_end     = 31225646;
my $gene_strand  = 1;
my $gene_name2   = "accC";
my $gene_start2  = 31225746;
my $gene_end2    = 31225946;
my $gene_strand2 = 1;

# create a gene/transcript/translation/operon set and then store it
my $gene        = Bio::EnsEMBL::Gene->new();
my $transcript  = Bio::EnsEMBL::Transcript->new();
my $exon        = Bio::EnsEMBL::Exon->new();
my $translation = Bio::EnsEMBL::Translation->new();
my $analysis    = $dba->get_AnalysisAdaptor->fetch_by_logic_name("Genscan");
ok( defined $analysis );
for my $feature ( $gene, $transcript, $exon ) {
	$feature->start($gene_start);
	$feature->end($gene_end);
	$feature->slice($slice);
	$feature->strand($gene_strand);
}
$gene->analysis($analysis);
$transcript->analysis($analysis);
$exon->phase(0);
$exon->end_phase(0);
$transcript->add_Exon($exon);
$translation->start_Exon($exon);
$translation->end_Exon($exon);
$translation->start(1);
$translation->end( ( $gene_start - $gene_end + 1 )/3 );
$translation->transcript($transcript);
$transcript->translation($translation);
$gene->add_Transcript($transcript);
$dba->get_GeneAdaptor()->store($gene);
ok( defined $gene->dbID() );
ok( defined $transcript->dbID() );
ok( defined $translation->dbID() );
ok( defined $exon->dbID() );
# create a second gene/transcript/translation/operon set and then store it
my $gene2        = Bio::EnsEMBL::Gene->new();
my $transcript2  = Bio::EnsEMBL::Transcript->new();
my $exon2        = Bio::EnsEMBL::Exon->new();
my $translation2 = Bio::EnsEMBL::Translation->new();
for my $feature ( $gene2, $transcript2, $exon2 ) {
	$feature->start($gene_start);
	$feature->end($gene_end);
	$feature->slice($slice);
	$feature->strand($gene_strand);
}
$gene2->analysis($analysis);
$transcript2->analysis($analysis);
$exon2->phase(0);
$exon2->end_phase(0);
$transcript2->add_Exon($exon2);
$translation2->start_Exon($exon2);
$translation2->end_Exon($exon2);
$translation2->start(1);
$translation2->end( ( $gene_start2 - $gene_end2 + 1 )/3 );
$translation2->transcript($transcript2);
$transcript2->translation($translation2);
$gene2->add_Transcript($transcript2);
$dba->get_GeneAdaptor()->store($gene2);
ok( defined $gene2->dbID() );
ok( defined $transcript2->dbID() );
ok( defined $translation2->dbID() );
ok( defined $exon2->dbID() );
# add the gene to the operon transcript
my $operon_adaptor = Bio::EnsEMBL::DBSQL::OperonAdaptor->new($dba);

# add the operon_transcript to the operon

# create an operon transcript
my $operon_transcript =
  Bio::EnsEMBL::OperonTranscript->new( -START  => $start,
									   -END    => $end,
									   -STRAND => $strand,
									   -SLICE  => $slice, -STABLE_ID=>"opt1" );
$operon_transcript->add_gene($gene);
$operon_transcript->add_gene($gene2);
$operon->add_OperonTranscript($operon_transcript);

my $operon_transcript2 =
  Bio::EnsEMBL::OperonTranscript->new( -START  => $start,
									   -END    => $gene_start2,
									   -STRAND => $strand,
									   -SLICE  => $slice, -STABLE_ID=>"opt2" );
$operon_transcript2->add_gene($gene);
$operon->add_OperonTranscript($operon_transcript2);

# store the lot
# store operon
$operon_adaptor->store($operon);
ok( defined $operon->dbID() );

# retrieve operon
my $operon2 = $operon_adaptor->fetch_by_dbID( $operon->dbID() );
ok( $operon2->dbID(),             $operon->dbID(),             "Operon ID" );
ok( $operon2->display_label(),    $operon->display_label(),    "Operon name" );
ok( $operon2->seq_region_start(), $operon->seq_region_start(), "Operon start" );
ok( $operon2->seq_region_end(),   $operon->seq_region_end(),   "Operon end" );
ok( $operon2->seq_region_strand(),
	$operon->seq_region_strand(),
	"Operon strand" );
my @operon_transcripts = @{ $operon2->get_all_OperonTranscripts() };
# check operon transcript
ok( scalar @operon_transcripts, 2, "Expected number of transcripts" );

my $operon_transcript_r = $operon_transcripts[0];
ok( $operon_transcript_r->dbID(), $operon_transcript->dbID(),
	"Operon transcript ID" );
ok( $operon_transcript_r->seq_region_start(),
	$operon_transcript->seq_region_start(),
	"Operon transcript start" );
ok( $operon_transcript_r->seq_region_end(),
	$operon_transcript->seq_region_end(),
	"Operon transcript end" );
ok( $operon_transcript_r->seq_region_strand(),
	$operon_transcript->seq_region_strand(),
	"Operon transcript strand" );

# check genes
my @ogenes = @{ $operon_transcript_r->get_all_Genes() };
ok( scalar @ogenes, 2, "Expected number of genes" );
my $gene_r = $ogenes[0];
ok( $gene_r->dbID(),              $gene->dbID(),              "Gene ID" );
ok( $gene_r->seq_region_start(),  $gene->seq_region_start(),  "Gene start" );
ok( $gene_r->seq_region_end(),    $gene->seq_region_end(),    "Gene end" );
ok( $gene_r->seq_region_strand(), $gene->seq_region_strand(), "Gene strand" );
$gene_r = $ogenes[1];
ok( $gene_r->dbID(),$gene2->dbID(), "Gene ID" );
ok( $gene_r->seq_region_start(),  $gene2->seq_region_start(),  "Gene start" );
ok( $gene_r->seq_region_end(),    $gene2->seq_region_end(),    "Gene end" );
ok( $gene_r->seq_region_strand(), $gene2->seq_region_strand(), "Gene strand" );

my $operon_transcript_r2 = $operon_transcripts[1];
ok( $operon_transcript_r2->dbID(),
	$operon_transcript2->dbID(),
	"Operon transcript ID" );
ok( $operon_transcript_r2->seq_region_start(),
	$operon_transcript2->seq_region_start(),
	"Operon transcript start" );
ok( $operon_transcript_r2->seq_region_end(),
	$operon_transcript2->seq_region_end(),
	"Operon transcript end" );
ok( $operon_transcript_r2->seq_region_strand(),
	$operon_transcript2->seq_region_strand(),
	"Operon transcript strand" );

# check genes
@ogenes = @{ $operon_transcript_r2->get_all_Genes() };
ok( scalar @ogenes, 1, "Expected number of genes" );
$gene_r = $ogenes[0];
ok( $gene_r->dbID(),              $gene->dbID(),              "Gene ID" );
ok( $gene_r->seq_region_start(),  $gene->seq_region_start(),  "Gene start" );
ok( $gene_r->seq_region_end(),    $gene->seq_region_end(),    "Gene end" );
ok( $gene_r->seq_region_strand(), $gene->seq_region_strand(), "Gene strand" );

$dba->get_OperonAdaptor()->remove($operon);
$dba->get_GeneAdaptor()->remove($gene);
$dba->get_GeneAdaptor()->remove($gene2);

