# script to convert 120 database to 130 database
# uses SchemaConvert
# see doc there


use SchemaConverter;
use DBI;

my $sourcedbh = DBI->connect("dbi:mysql:host=ecs2d;database=homo_sapiens_core_6_28", "ensro");
my $targetdbh = DBI->connect("dbi:mysql:host=ecs1c.sanger.ac.uk;database=arne_main_28", "ensadmin", "ensembl");

my $sc = SchemaConverter->new( $sourcedbh, $targetdbh );
$sc->tmp_dir( "/acari/work1/stabenau/db" );


$sc->custom_select( "assembly", "select s.fpcctg_name, c.chromosome_id, s.raw_id, s.chr_start, s.chr_end, s.fpcctg_start, s.fpcctg_end, s.fpcctg_ori, s.raw_start, s.raw_end, s.raw_ori, s.type from chromosome c, static_golden_path s where s.chr_name = c.name" ); 

$sc->table_rename( "analysisprocess", "analysis" );
$sc->column_rename( "analysis", "analysisId", "analysis_id" );

$sc->column_rename(  "contig", "internal_id", "contig_id" );
$sc->column_rename(  "contig", "id", "name" );
$sc->column_rename(  "contig", "clone", "clone_id" );
$sc->column_rename(  "contig", "dna", "dna_id" );
$sc->column_rename(  "contig", "chromosomeId", "chromosome_id" );
$sc->column_rename(  "contig", "international_id", "international_name" );

$sc->column_rename(  "clone", "internal_id", "clone_id" );
$sc->column_rename(  "clone", "id", "name" );
$sc->column_rename(  "clone", "embl_id", "embl_acc" );

$sc->custom_select(  "map_density", "select c.chromosome_id, m.chr_start, m.chr_end, m.type, m.value from chromosome c, map_density m where m.chr_name = c.name" ); 
$sc->custom_select(  "karyotype", "select c.chromosome_id, k.chr_start, k.chr_end, k.band, k.stain from chromosome c, karyotype k where k.chr_name = c.name" ); 


$sc->column_rename(  "dna", "id", "dna_id" );
$sc->big_table( "dna" );

$sc->table_skip(  "exon_feature" );
$sc->table_skip(  "simple_feature" );
$sc->table_skip(  "dna_align_feature" );
$sc->table_skip(  "protein_align_feature" );
$sc->table_skip(  "repeat_feature" );
$sc->table_skip(  "repeat_consensus" );
$sc->table_skip(  "assembly_locations" );
$sc->table_skip( "landmark_marker" );
$sc->table_skip( "prediction_transcript" );
$sc->table_skip( "supporting_feature" );

$sc->table_rename(  "objectXref", "object_xref" );
$sc->column_rename(   "object_xref", "objectxrefId","object_xref_id" );
$sc->column_rename(   "object_xref", "xrefId", "xref_id" );

$sc->table_rename(  "identityXref", "identity_xref" );
$sc->column_rename(   "identity_xref", "objectxrefId","object_xref_id" );

$sc->table_rename(  "Xref", "xref" );
$sc->column_rename(  "xref", "xrefId", "xref_id" );
$sc->column_rename(  "xref", "externalDBId", "external_db_id" );
$sc->column_rename(  "xref", "dbprimary_id", "dbprimary_acc" );
$sc->column_rename(  "xref", "display_id", "display_label" );

$sc->table_rename(  "externalSynonym", "external_synonym" );
$sc->column_rename(  "external_synonym", "xrefId", "xref_id" );

$sc->table_rename(  "externalDB", "external_db" );
$sc->column_rename(  "external_db", "externalDBId", "external_db_id" );

$sc->column_skip(  "supporting_feature", "contig_id" );

$sc->column_rename(  "protein_feature", "id", "protein_feature_id" );
$sc->column_rename(  "protein_feature", "translation", "translation_id" );
$sc->column_rename(  "protein_feature", "analysis", "analysis_id" );
$sc->column_rename(  "protein_feature", "hstart", "hit_start" );
$sc->column_rename(  "protein_feature", "hend", "hit_end" );
$sc->column_rename(  "protein_feature", "hid", "hit_id" );
$sc->column_rename(  "protein_feature", "perc_id", "perc_ident" );

$sc->column_rename(  "exon", "seq_start", "contig_start" );
$sc->column_rename(  "exon", "seq_end", "contig_end" );
$sc->column_rename(  "exon", "strand", "contig_strand" );

$sc->column_rename(  "gene", "analysisId", "analysis_id" );
$sc->column_skip( "transcript","exon_count" );
$sc->column_skip( "gene", "transcript_count" );

$sc->column_rename(  "supporting_feature", "hid", "hit_id" );
$sc->column_rename(  "supporting_feature", "hstart", "hit_start" );
$sc->column_rename(  "supporting_feature", "hend", "hit_end" );
$sc->column_rename(  "supporting_feature", "seq_start", "contig_start" );
$sc->column_rename(  "supporting_feature", "seq_end", "contig_end" );
$sc->column_rename(  "supporting_feature", "hstrand", "hit_strand" );
$sc->column_rename(  "supporting_feature", "analysis", "analysis_id" );

$sc->clear_target();
$sc->transfer();



