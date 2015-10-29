-- 
-- Created by SQL::Translator::Producer::SQLite
-- Created on Thu Oct 29 15:38:20 2015
-- 

BEGIN TRANSACTION;

--
-- Table: alt_allele
--
CREATE TABLE alt_allele (
  alt_allele_id INTEGER PRIMARY KEY NOT NULL,
  alt_allele_group_id int(10) NOT NULL,
  gene_id int(10) NOT NULL
);

CREATE INDEX gene_id ON alt_allele (gene_id, alt_allele_group_id);

CREATE UNIQUE INDEX gene_idx ON alt_allele (gene_id);

--
-- Table: alt_allele_attrib
--
CREATE TABLE alt_allele_attrib (
  alt_allele_id int(10) DEFAULT NULL,
  attrib enum(35) DEFAULT NULL
);

CREATE INDEX aa_idx ON alt_allele_attrib (alt_allele_id, attrib);

--
-- Table: alt_allele_group
--
CREATE TABLE alt_allele_group (
  alt_allele_group_id INTEGER PRIMARY KEY NOT NULL
);

--
-- Table: analysis
--
CREATE TABLE analysis (
  analysis_id INTEGER PRIMARY KEY NOT NULL,
  created datetime NOT NULL DEFAULT '0000-00-00 00:00:00',
  logic_name varchar(128) NOT NULL,
  db varchar(120) DEFAULT NULL,
  db_version varchar(40) DEFAULT NULL,
  db_file varchar(120) DEFAULT NULL,
  program varchar(80) DEFAULT NULL,
  program_version varchar(40) DEFAULT NULL,
  program_file varchar(80) DEFAULT NULL,
  parameters text,
  module varchar(80) DEFAULT NULL,
  module_version varchar(40) DEFAULT NULL,
  gff_source varchar(40) DEFAULT NULL,
  gff_feature varchar(40) DEFAULT NULL
);

CREATE UNIQUE INDEX logic_name_idx ON analysis (logic_name);

--
-- Table: analysis_description
--
CREATE TABLE analysis_description (
  analysis_id smallint(5) NOT NULL,
  description text,
  display_label varchar(255) NOT NULL,
  displayable tinyint(1) NOT NULL DEFAULT 1,
  web_data text
);

CREATE UNIQUE INDEX analysis_idx ON analysis_description (analysis_id);

--
-- Table: assembly
--
CREATE TABLE assembly (
  asm_seq_region_id int(10) NOT NULL,
  cmp_seq_region_id int(10) NOT NULL,
  asm_start int(10) NOT NULL,
  asm_end int(10) NOT NULL,
  cmp_start int(10) NOT NULL,
  cmp_end int(10) NOT NULL,
  ori tinyint(4) NOT NULL
);

CREATE INDEX cmp_seq_region_idx ON assembly (cmp_seq_region_id);

CREATE INDEX asm_seq_region_idx ON assembly (asm_seq_region_id, asm_start);

CREATE UNIQUE INDEX all_idx ON assembly (asm_seq_region_id, cmp_seq_region_id, asm_start, asm_end, cmp_start, cmp_end, ori);

--
-- Table: assembly_exception
--
CREATE TABLE assembly_exception (
  assembly_exception_id INTEGER PRIMARY KEY NOT NULL,
  seq_region_id int(10) NOT NULL,
  seq_region_start int(10) NOT NULL,
  seq_region_end int(10) NOT NULL,
  exc_type enum(11) NOT NULL,
  exc_seq_region_id int(10) NOT NULL,
  exc_seq_region_start int(10) NOT NULL,
  exc_seq_region_end int(10) NOT NULL,
  ori int(11) NOT NULL
);

CREATE INDEX sr_idx ON assembly_exception (seq_region_id, seq_region_start);

CREATE INDEX ex_idx ON assembly_exception (exc_seq_region_id, exc_seq_region_start);

--
-- Table: associated_group
--
CREATE TABLE associated_group (
  associated_group_id INTEGER PRIMARY KEY NOT NULL,
  description varchar(128) DEFAULT NULL
);

--
-- Table: associated_xref
--
CREATE TABLE associated_xref (
  associated_xref_id INTEGER PRIMARY KEY NOT NULL,
  object_xref_id int(10) NOT NULL DEFAULT 0,
  xref_id int(10) NOT NULL DEFAULT 0,
  source_xref_id int(10) DEFAULT NULL,
  condition_type varchar(128) DEFAULT NULL,
  associated_group_id int(10) DEFAULT NULL,
  rank int(10) DEFAULT 0
);

CREATE INDEX associated_source_idx ON associated_xref (source_xref_id);

CREATE INDEX associated_object_idx ON associated_xref (object_xref_id);

CREATE INDEX associated_idx ON associated_xref (xref_id);

CREATE INDEX associated_group_idx ON associated_xref (associated_group_id);

CREATE UNIQUE INDEX object_associated_source_type_idx ON associated_xref (object_xref_id, xref_id, source_xref_id, condition_type, associated_group_id);

--
-- Table: attrib_type
--
CREATE TABLE attrib_type (
  attrib_type_id INTEGER PRIMARY KEY NOT NULL,
  code varchar(20) NOT NULL DEFAULT '',
  name varchar(255) NOT NULL DEFAULT '',
  description text
);

CREATE UNIQUE INDEX code_idx ON attrib_type (code);

--
-- Table: coord_system
--
CREATE TABLE coord_system (
  coord_system_id INTEGER PRIMARY KEY NOT NULL,
  species_id int(10) NOT NULL DEFAULT 1,
  name varchar(40) NOT NULL,
  version varchar(255) DEFAULT NULL,
  rank int(11) NOT NULL,
  attrib varchar(15) DEFAULT NULL
);

CREATE INDEX species_idx ON coord_system (species_id);

CREATE UNIQUE INDEX rank_idx ON coord_system (rank, species_id);

CREATE UNIQUE INDEX name_idx ON coord_system (name, version, species_id);

--
-- Table: data_file
--
CREATE TABLE data_file (
  data_file_id INTEGER PRIMARY KEY NOT NULL,
  coord_system_id int(10) NOT NULL,
  analysis_id smallint(5) NOT NULL,
  name varchar(100) NOT NULL,
  version_lock tinyint(1) NOT NULL DEFAULT 0,
  absolute tinyint(1) NOT NULL DEFAULT 0,
  url text,
  file_type enum(6) DEFAULT NULL
);

CREATE INDEX df_name_idx ON data_file (name);

CREATE INDEX df_analysis_idx ON data_file (analysis_id);

CREATE UNIQUE INDEX df_unq_idx ON data_file (coord_system_id, analysis_id, name, file_type);

--
-- Table: density_feature
--
CREATE TABLE density_feature (
  density_feature_id INTEGER PRIMARY KEY NOT NULL,
  density_type_id int(10) NOT NULL,
  seq_region_id int(10) NOT NULL,
  seq_region_start int(10) NOT NULL,
  seq_region_end int(10) NOT NULL,
  density_value float(8,2) NOT NULL
);

CREATE INDEX seq_region_idx ON density_feature (density_type_id, seq_region_id, seq_region_start);

CREATE INDEX seq_region_id_idx ON density_feature (seq_region_id);

--
-- Table: density_type
--
CREATE TABLE density_type (
  density_type_id INTEGER PRIMARY KEY NOT NULL,
  analysis_id smallint(5) NOT NULL,
  block_size int(11) NOT NULL,
  region_features int(11) NOT NULL,
  value_type enum(5) NOT NULL
);

CREATE UNIQUE INDEX analysis_idx02 ON density_type (analysis_id, block_size, region_features);

--
-- Table: dependent_xref
--
CREATE TABLE dependent_xref (
  object_xref_id INTEGER PRIMARY KEY NOT NULL,
  master_xref_id int(10) NOT NULL,
  dependent_xref_id int(10) NOT NULL
);

CREATE INDEX dependent ON dependent_xref (dependent_xref_id);

CREATE INDEX master_idx ON dependent_xref (master_xref_id);

--
-- Table: ditag
--
CREATE TABLE ditag (
  ditag_id INTEGER PRIMARY KEY NOT NULL,
  name varchar(30) NOT NULL DEFAULT '',
  type varchar(30) NOT NULL DEFAULT '',
  tag_count smallint(6) NOT NULL DEFAULT 1,
  sequence tinytext(255) NOT NULL
);

--
-- Table: ditag_feature
--
CREATE TABLE ditag_feature (
  ditag_feature_id INTEGER PRIMARY KEY NOT NULL,
  ditag_id int(10) NOT NULL DEFAULT 0,
  ditag_pair_id int(10) NOT NULL DEFAULT 0,
  seq_region_id int(10) NOT NULL DEFAULT 0,
  seq_region_start int(10) NOT NULL DEFAULT 0,
  seq_region_end int(10) NOT NULL DEFAULT 0,
  seq_region_strand tinyint(1) NOT NULL DEFAULT 0,
  analysis_id smallint(5) NOT NULL DEFAULT 0,
  hit_start int(10) NOT NULL DEFAULT 0,
  hit_end int(10) NOT NULL DEFAULT 0,
  hit_strand tinyint(1) NOT NULL DEFAULT 0,
  cigar_line tinytext(255) NOT NULL,
  ditag_side enum(1) NOT NULL
);

CREATE INDEX seq_region_idx02 ON ditag_feature (seq_region_id, seq_region_start, seq_region_end);

CREATE INDEX ditag_idx ON ditag_feature (ditag_id);

CREATE INDEX ditag_pair_idx ON ditag_feature (ditag_pair_id);

--
-- Table: dna
--
CREATE TABLE dna (
  seq_region_id INTEGER PRIMARY KEY NOT NULL,
  sequence longtext(4294967295) NOT NULL
);

--
-- Table: dna_align_feature
--
CREATE TABLE dna_align_feature (
  dna_align_feature_id INTEGER PRIMARY KEY NOT NULL,
  seq_region_id int(10) NOT NULL,
  seq_region_start int(10) NOT NULL,
  seq_region_end int(10) NOT NULL,
  seq_region_strand tinyint(1) NOT NULL,
  hit_start int(11) NOT NULL,
  hit_end int(11) NOT NULL,
  hit_strand tinyint(1) NOT NULL,
  hit_name varchar(40) NOT NULL,
  analysis_id smallint(5) NOT NULL,
  score double(8,2) DEFAULT NULL,
  evalue double(8,2) DEFAULT NULL,
  perc_ident float(8,2) DEFAULT NULL,
  cigar_line text,
  external_db_id int(10) DEFAULT NULL,
  hcoverage double(8,2) DEFAULT NULL,
  external_data text
);

CREATE INDEX seq_region_idx03 ON dna_align_feature (seq_region_id, analysis_id, seq_region_start, score);

CREATE INDEX seq_region_idx_2 ON dna_align_feature (seq_region_id, seq_region_start);

CREATE INDEX hit_idx ON dna_align_feature (hit_name);

CREATE INDEX analysis_idx03 ON dna_align_feature (analysis_id);

CREATE INDEX external_db_idx ON dna_align_feature (external_db_id);

--
-- Table: exon
--
CREATE TABLE exon (
  exon_id INTEGER PRIMARY KEY NOT NULL,
  seq_region_id int(10) NOT NULL,
  seq_region_start int(10) NOT NULL,
  seq_region_end int(10) NOT NULL,
  seq_region_strand tinyint(2) NOT NULL,
  phase tinyint(2) NOT NULL,
  end_phase tinyint(2) NOT NULL,
  is_current tinyint(1) NOT NULL DEFAULT 1,
  is_constitutive tinyint(1) NOT NULL DEFAULT 0,
  stable_id varchar(128) DEFAULT NULL,
  version smallint(5) NOT NULL DEFAULT 1,
  created_date datetime NOT NULL DEFAULT '0000-00-00 00:00:00',
  modified_date datetime NOT NULL DEFAULT '0000-00-00 00:00:00'
);

CREATE INDEX seq_region_idx04 ON exon (seq_region_id, seq_region_start);

CREATE INDEX stable_id_idx ON exon (stable_id, version);

--
-- Table: exon_transcript
--
CREATE TABLE exon_transcript (
  exon_id int(10) NOT NULL,
  transcript_id int(10) NOT NULL,
  rank int(10) NOT NULL,
  PRIMARY KEY (exon_id, transcript_id, rank)
);

CREATE INDEX transcript ON exon_transcript (transcript_id);

CREATE INDEX exon02 ON exon_transcript (exon_id);

--
-- Table: external_db
--
CREATE TABLE external_db (
  external_db_id INTEGER PRIMARY KEY NOT NULL,
  db_name varchar(100) NOT NULL,
  db_release varchar(255) DEFAULT NULL,
  status enum(9) NOT NULL,
  dbprimary_acc_linkable tinyint(1) NOT NULL DEFAULT 1,
  display_label_linkable tinyint(1) NOT NULL DEFAULT 0,
  priority int(11) NOT NULL,
  db_display_name varchar(255) DEFAULT NULL,
  type enum(18) DEFAULT NULL,
  secondary_db_name varchar(255) DEFAULT NULL,
  secondary_db_table varchar(255) DEFAULT NULL,
  description text
);

--
-- Table: external_synonym
--
CREATE TABLE external_synonym (
  xref_id int(10) NOT NULL,
  synonym varchar(100) NOT NULL,
  PRIMARY KEY (xref_id, synonym)
);

CREATE INDEX name_index ON external_synonym (synonym);

--
-- Table: gene
--
CREATE TABLE gene (
  gene_id INTEGER PRIMARY KEY NOT NULL,
  biotype varchar(40) NOT NULL,
  analysis_id smallint(5) NOT NULL,
  seq_region_id int(10) NOT NULL,
  seq_region_start int(10) NOT NULL,
  seq_region_end int(10) NOT NULL,
  seq_region_strand tinyint(2) NOT NULL,
  display_xref_id int(10) DEFAULT NULL,
  source varchar(40) NOT NULL,
  status enum(19) DEFAULT NULL,
  description text,
  is_current tinyint(1) NOT NULL DEFAULT 1,
  canonical_transcript_id int(10) NOT NULL,
  stable_id varchar(128) DEFAULT NULL,
  version smallint(5) NOT NULL DEFAULT 1,
  created_date datetime NOT NULL DEFAULT '0000-00-00 00:00:00',
  modified_date datetime NOT NULL DEFAULT '0000-00-00 00:00:00'
);

CREATE INDEX seq_region_idx05 ON gene (seq_region_id, seq_region_start);

CREATE INDEX xref_id_index ON gene (display_xref_id);

CREATE INDEX analysis_idx04 ON gene (analysis_id);

CREATE INDEX stable_id_idx02 ON gene (stable_id, version);

CREATE INDEX canonical_transcript_id_idx ON gene (canonical_transcript_id);

--
-- Table: gene_archive
--
CREATE TABLE gene_archive (
  gene_stable_id varchar(128) NOT NULL,
  gene_version smallint(6) NOT NULL DEFAULT 1,
  transcript_stable_id varchar(128) NOT NULL,
  transcript_version smallint(6) NOT NULL DEFAULT 1,
  translation_stable_id varchar(128) DEFAULT NULL,
  translation_version smallint(6) NOT NULL DEFAULT 1,
  peptide_archive_id int(10) DEFAULT NULL,
  mapping_session_id int(10) NOT NULL
);

CREATE INDEX peptide_archive_id_idx ON gene_archive (peptide_archive_id);

CREATE INDEX gene_idx02 ON gene_archive (gene_stable_id, gene_version);

CREATE INDEX transcript_idx ON gene_archive (transcript_stable_id, transcript_version);

CREATE INDEX translation_idx ON gene_archive (translation_stable_id, translation_version);

--
-- Table: gene_attrib
--
CREATE TABLE gene_attrib (
  gene_id int(10) NOT NULL DEFAULT 0,
  attrib_type_id smallint(5) NOT NULL DEFAULT 0,
  value text NOT NULL
);

CREATE INDEX type_val_idx ON gene_attrib (attrib_type_id, value);

CREATE INDEX val_only_idx ON gene_attrib (value);

CREATE INDEX gene_idx03 ON gene_attrib (gene_id);

--
-- Table: genome_statistics
--
CREATE TABLE genome_statistics (
  genome_statistics_id INTEGER PRIMARY KEY NOT NULL,
  statistic varchar(128) NOT NULL,
  value bigint(11) NOT NULL DEFAULT 0,
  species_id int(10) DEFAULT 1,
  attrib_type_id int(10) DEFAULT NULL,
  timestamp datetime NOT NULL DEFAULT '0000-00-00 00:00:00'
);

CREATE INDEX stats_idx ON genome_statistics (statistic, attrib_type_id, species_id);

CREATE UNIQUE INDEX stats_uniq ON genome_statistics (statistic, attrib_type_id, species_id);

--
-- Table: identity_xref
--
CREATE TABLE identity_xref (
  object_xref_id INTEGER PRIMARY KEY NOT NULL,
  xref_identity int(5) DEFAULT NULL,
  ensembl_identity int(5) DEFAULT NULL,
  xref_start int(11) DEFAULT NULL,
  xref_end int(11) DEFAULT NULL,
  ensembl_start int(11) DEFAULT NULL,
  ensembl_end int(11) DEFAULT NULL,
  cigar_line text,
  score double(8,2) DEFAULT NULL,
  evalue double(8,2) DEFAULT NULL
);

--
-- Table: interpro
--
CREATE TABLE interpro (
  interpro_ac varchar(40) NOT NULL,
  id varchar(40) NOT NULL
);

CREATE INDEX id_idx ON interpro (id);

CREATE UNIQUE INDEX accession_idx ON interpro (interpro_ac, id);

--
-- Table: intron_supporting_evidence
--
CREATE TABLE intron_supporting_evidence (
  intron_supporting_evidence_id INTEGER PRIMARY KEY NOT NULL,
  analysis_id smallint(5) NOT NULL,
  seq_region_id int(10) NOT NULL,
  seq_region_start int(10) NOT NULL,
  seq_region_end int(10) NOT NULL,
  seq_region_strand tinyint(2) NOT NULL,
  hit_name varchar(100) NOT NULL,
  score decimal(10,3) DEFAULT NULL,
  score_type enum(5) DEFAULT 'NONE',
  is_splice_canonical tinyint(1) NOT NULL DEFAULT 0
);

CREATE INDEX seq_region_idx06 ON intron_supporting_evidence (seq_region_id, seq_region_start);

CREATE UNIQUE INDEX analysis_id ON intron_supporting_evidence (analysis_id, seq_region_id, seq_region_start, seq_region_end, seq_region_strand, hit_name);

--
-- Table: karyotype
--
CREATE TABLE karyotype (
  karyotype_id INTEGER PRIMARY KEY NOT NULL,
  seq_region_id int(10) NOT NULL,
  seq_region_start int(10) NOT NULL,
  seq_region_end int(10) NOT NULL,
  band varchar(40) DEFAULT NULL,
  stain varchar(40) DEFAULT NULL
);

CREATE INDEX region_band_idx ON karyotype (seq_region_id, band);

--
-- Table: map
--
CREATE TABLE map (
  map_id INTEGER PRIMARY KEY NOT NULL,
  map_name varchar(30) NOT NULL
);

--
-- Table: mapping_session
--
CREATE TABLE mapping_session (
  mapping_session_id INTEGER PRIMARY KEY NOT NULL,
  old_db_name varchar(80) NOT NULL DEFAULT '',
  new_db_name varchar(80) NOT NULL DEFAULT '',
  old_release varchar(5) NOT NULL DEFAULT '',
  new_release varchar(5) NOT NULL DEFAULT '',
  old_assembly varchar(20) NOT NULL DEFAULT '',
  new_assembly varchar(20) NOT NULL DEFAULT '',
  created datetime NOT NULL
);

--
-- Table: mapping_set
--
CREATE TABLE mapping_set (
  mapping_set_id INTEGER PRIMARY KEY NOT NULL,
  internal_schema_build varchar(20) NOT NULL,
  external_schema_build varchar(20) NOT NULL
);

CREATE UNIQUE INDEX mapping_idx ON mapping_set (internal_schema_build, external_schema_build);

--
-- Table: marker
--
CREATE TABLE marker (
  marker_id INTEGER PRIMARY KEY NOT NULL,
  display_marker_synonym_id int(10) DEFAULT NULL,
  left_primer varchar(100) NOT NULL,
  right_primer varchar(100) NOT NULL,
  min_primer_dist int(10) NOT NULL,
  max_primer_dist int(10) NOT NULL,
  priority int(11) DEFAULT NULL,
  type enum(14) DEFAULT NULL
);

CREATE INDEX marker_idx ON marker (marker_id, priority);

CREATE INDEX display_idx ON marker (display_marker_synonym_id);

--
-- Table: marker_feature
--
CREATE TABLE marker_feature (
  marker_feature_id INTEGER PRIMARY KEY NOT NULL,
  marker_id int(10) NOT NULL,
  seq_region_id int(10) NOT NULL,
  seq_region_start int(10) NOT NULL,
  seq_region_end int(10) NOT NULL,
  analysis_id smallint(5) NOT NULL,
  map_weight int(10) DEFAULT NULL
);

CREATE INDEX seq_region_idx07 ON marker_feature (seq_region_id, seq_region_start);

CREATE INDEX analysis_idx05 ON marker_feature (analysis_id);

--
-- Table: marker_map_location
--
CREATE TABLE marker_map_location (
  marker_id int(10) NOT NULL,
  map_id int(10) NOT NULL,
  chromosome_name varchar(15) NOT NULL,
  marker_synonym_id int(10) NOT NULL,
  position varchar(15) NOT NULL,
  lod_score double(8,2) DEFAULT NULL,
  PRIMARY KEY (marker_id, map_id)
);

CREATE INDEX map_idx ON marker_map_location (map_id, chromosome_name, position);

--
-- Table: marker_synonym
--
CREATE TABLE marker_synonym (
  marker_synonym_id INTEGER PRIMARY KEY NOT NULL,
  marker_id int(10) NOT NULL,
  source varchar(20) DEFAULT NULL,
  name varchar(50) DEFAULT NULL
);

CREATE INDEX marker_synonym_idx ON marker_synonym (marker_synonym_id, name);

CREATE INDEX marker_idx02 ON marker_synonym (marker_id);

--
-- Table: meta
--
CREATE TABLE meta (
  meta_id INTEGER PRIMARY KEY NOT NULL,
  species_id int(10) DEFAULT 1,
  meta_key varchar(40) NOT NULL,
  meta_value varchar(255) DEFAULT NULL
);

CREATE INDEX species_value_idx ON meta (species_id, meta_value);

CREATE UNIQUE INDEX species_key_value_idx ON meta (species_id, meta_key, meta_value);

--
-- Table: meta_coord
--
CREATE TABLE meta_coord (
  table_name varchar(40) NOT NULL,
  coord_system_id int(10) NOT NULL,
  max_length int(11) DEFAULT NULL
);

CREATE UNIQUE INDEX cs_table_name_idx ON meta_coord (coord_system_id, table_name);

--
-- Table: misc_attrib
--
CREATE TABLE misc_attrib (
  misc_feature_id int(10) NOT NULL DEFAULT 0,
  attrib_type_id smallint(5) NOT NULL DEFAULT 0,
  value text NOT NULL
);

CREATE INDEX type_val_idx02 ON misc_attrib (attrib_type_id, value);

CREATE INDEX val_only_idx02 ON misc_attrib (value);

CREATE INDEX misc_feature_idx ON misc_attrib (misc_feature_id);

--
-- Table: misc_feature
--
CREATE TABLE misc_feature (
  misc_feature_id INTEGER PRIMARY KEY NOT NULL,
  seq_region_id int(10) NOT NULL DEFAULT 0,
  seq_region_start int(10) NOT NULL DEFAULT 0,
  seq_region_end int(10) NOT NULL DEFAULT 0,
  seq_region_strand tinyint(4) NOT NULL DEFAULT 0
);

CREATE INDEX seq_region_idx08 ON misc_feature (seq_region_id, seq_region_start);

--
-- Table: misc_feature_misc_set
--
CREATE TABLE misc_feature_misc_set (
  misc_feature_id int(10) NOT NULL DEFAULT 0,
  misc_set_id smallint(5) NOT NULL DEFAULT 0,
  PRIMARY KEY (misc_feature_id, misc_set_id)
);

CREATE INDEX reverse_idx ON misc_feature_misc_set (misc_set_id, misc_feature_id);

--
-- Table: misc_set
--
CREATE TABLE misc_set (
  misc_set_id INTEGER PRIMARY KEY NOT NULL,
  code varchar(25) NOT NULL DEFAULT '',
  name varchar(255) NOT NULL DEFAULT '',
  description text NOT NULL,
  max_length int(10) NOT NULL
);

CREATE UNIQUE INDEX code_idx02 ON misc_set (code);

--
-- Table: object_xref
--
CREATE TABLE object_xref (
  object_xref_id INTEGER PRIMARY KEY NOT NULL,
  ensembl_id int(10) NOT NULL,
  ensembl_object_type enum(16) NOT NULL,
  xref_id int(10) NOT NULL,
  linkage_annotation varchar(255) DEFAULT NULL,
  analysis_id smallint(5) NOT NULL DEFAULT 0
);

CREATE INDEX analysis_idx06 ON object_xref (analysis_id);

CREATE INDEX ensembl_idx ON object_xref (ensembl_object_type, ensembl_id);

CREATE UNIQUE INDEX xref_idx ON object_xref (xref_id, ensembl_object_type, ensembl_id, analysis_id);

--
-- Table: ontology_xref
--
CREATE TABLE ontology_xref (
  object_xref_id int(10) NOT NULL DEFAULT 0,
  source_xref_id int(10) DEFAULT NULL,
  linkage_type varchar(3) DEFAULT NULL
);

CREATE INDEX source_idx ON ontology_xref (source_xref_id);

CREATE UNIQUE INDEX object_source_type_idx ON ontology_xref (object_xref_id, source_xref_id, linkage_type);

--
-- Table: operon
--
CREATE TABLE operon (
  operon_id INTEGER PRIMARY KEY NOT NULL,
  seq_region_id int(10) NOT NULL,
  seq_region_start int(10) NOT NULL,
  seq_region_end int(10) NOT NULL,
  seq_region_strand tinyint(2) NOT NULL,
  display_label varchar(255) DEFAULT NULL,
  analysis_id smallint(5) NOT NULL,
  stable_id varchar(128) DEFAULT NULL,
  version smallint(5) NOT NULL DEFAULT 1,
  created_date datetime NOT NULL DEFAULT '0000-00-00 00:00:00',
  modified_date datetime NOT NULL DEFAULT '0000-00-00 00:00:00'
);

CREATE INDEX seq_region_idx09 ON operon (seq_region_id, seq_region_start);

CREATE INDEX name_idx02 ON operon (display_label);

CREATE INDEX stable_id_idx03 ON operon (stable_id, version);

--
-- Table: operon_transcript
--
CREATE TABLE operon_transcript (
  operon_transcript_id INTEGER PRIMARY KEY NOT NULL,
  seq_region_id int(10) NOT NULL,
  seq_region_start int(10) NOT NULL,
  seq_region_end int(10) NOT NULL,
  seq_region_strand tinyint(2) NOT NULL,
  operon_id int(10) NOT NULL,
  display_label varchar(255) DEFAULT NULL,
  analysis_id smallint(5) NOT NULL,
  stable_id varchar(128) DEFAULT NULL,
  version smallint(5) NOT NULL DEFAULT 1,
  created_date datetime NOT NULL DEFAULT '0000-00-00 00:00:00',
  modified_date datetime NOT NULL DEFAULT '0000-00-00 00:00:00'
);

CREATE INDEX operon_idx ON operon_transcript (operon_id);

CREATE INDEX seq_region_idx10 ON operon_transcript (seq_region_id, seq_region_start);

CREATE INDEX stable_id_idx04 ON operon_transcript (stable_id, version);

--
-- Table: operon_transcript_gene
--
CREATE TABLE operon_transcript_gene (
  operon_transcript_id int(10) DEFAULT NULL,
  gene_id int(10) DEFAULT NULL
);

CREATE INDEX operon_transcript_gene_idx ON operon_transcript_gene (operon_transcript_id, gene_id);

--
-- Table: peptide_archive
--
CREATE TABLE peptide_archive (
  peptide_archive_id INTEGER PRIMARY KEY NOT NULL,
  md5_checksum varchar(32) DEFAULT NULL,
  peptide_seq mediumtext(16777215) NOT NULL
);

CREATE INDEX checksum ON peptide_archive (md5_checksum);

--
-- Table: prediction_exon
--
CREATE TABLE prediction_exon (
  prediction_exon_id INTEGER PRIMARY KEY NOT NULL,
  prediction_transcript_id int(10) NOT NULL,
  exon_rank smallint(5) NOT NULL,
  seq_region_id int(10) NOT NULL,
  seq_region_start int(10) NOT NULL,
  seq_region_end int(10) NOT NULL,
  seq_region_strand tinyint(4) NOT NULL,
  start_phase tinyint(4) NOT NULL,
  score double(8,2) DEFAULT NULL,
  p_value double(8,2) DEFAULT NULL
);

CREATE INDEX transcript_idx02 ON prediction_exon (prediction_transcript_id);

CREATE INDEX seq_region_idx11 ON prediction_exon (seq_region_id, seq_region_start);

--
-- Table: prediction_transcript
--
CREATE TABLE prediction_transcript (
  prediction_transcript_id INTEGER PRIMARY KEY NOT NULL,
  seq_region_id int(10) NOT NULL,
  seq_region_start int(10) NOT NULL,
  seq_region_end int(10) NOT NULL,
  seq_region_strand tinyint(4) NOT NULL,
  analysis_id smallint(5) NOT NULL,
  display_label varchar(255) DEFAULT NULL
);

CREATE INDEX analysis_idx07 ON prediction_transcript (analysis_id);

CREATE INDEX seq_region_idx12 ON prediction_transcript (seq_region_id, seq_region_start);

--
-- Table: protein_align_feature
--
CREATE TABLE protein_align_feature (
  protein_align_feature_id INTEGER PRIMARY KEY NOT NULL,
  seq_region_id int(10) NOT NULL,
  seq_region_start int(10) NOT NULL,
  seq_region_end int(10) NOT NULL,
  seq_region_strand tinyint(1) NOT NULL DEFAULT 1,
  hit_start int(10) NOT NULL,
  hit_end int(10) NOT NULL,
  hit_name varchar(40) NOT NULL,
  analysis_id smallint(5) NOT NULL,
  score double(8,2) DEFAULT NULL,
  evalue double(8,2) DEFAULT NULL,
  perc_ident float(8,2) DEFAULT NULL,
  cigar_line text,
  external_db_id int(10) DEFAULT NULL,
  hcoverage double(8,2) DEFAULT NULL
);

CREATE INDEX seq_region_idx13 ON protein_align_feature (seq_region_id, analysis_id, seq_region_start, score);

CREATE INDEX seq_region_idx_202 ON protein_align_feature (seq_region_id, seq_region_start);

CREATE INDEX hit_idx02 ON protein_align_feature (hit_name);

CREATE INDEX analysis_idx08 ON protein_align_feature (analysis_id);

CREATE INDEX external_db_idx02 ON protein_align_feature (external_db_id);

--
-- Table: protein_feature
--
CREATE TABLE protein_feature (
  protein_feature_id INTEGER PRIMARY KEY NOT NULL,
  translation_id int(10) NOT NULL,
  seq_start int(10) NOT NULL,
  seq_end int(10) NOT NULL,
  hit_start int(10) NOT NULL,
  hit_end int(10) NOT NULL,
  hit_name varchar(40) NOT NULL,
  analysis_id smallint(5) NOT NULL,
  score double(8,2) DEFAULT NULL,
  evalue double(8,2) DEFAULT NULL,
  perc_ident float(8,2) DEFAULT NULL,
  external_data text,
  hit_description text
);

CREATE INDEX hitname_idx ON protein_feature (hit_name);

CREATE INDEX analysis_idx09 ON protein_feature (analysis_id);

CREATE INDEX translation_idx02 ON protein_feature (translation_id);

--
-- Table: repeat_consensus
--
CREATE TABLE repeat_consensus (
  repeat_consensus_id INTEGER PRIMARY KEY NOT NULL,
  repeat_name varchar(255) NOT NULL,
  repeat_class varchar(100) NOT NULL,
  repeat_type varchar(40) NOT NULL,
  repeat_consensus text
);

CREATE INDEX name ON repeat_consensus (repeat_name);

CREATE INDEX class ON repeat_consensus (repeat_class);

CREATE INDEX consensus ON repeat_consensus (repeat_consensus);

CREATE INDEX type ON repeat_consensus (repeat_type);

--
-- Table: repeat_feature
--
CREATE TABLE repeat_feature (
  repeat_feature_id INTEGER PRIMARY KEY NOT NULL,
  seq_region_id int(10) NOT NULL,
  seq_region_start int(10) NOT NULL,
  seq_region_end int(10) NOT NULL,
  seq_region_strand tinyint(1) NOT NULL DEFAULT 1,
  repeat_start int(10) NOT NULL,
  repeat_end int(10) NOT NULL,
  repeat_consensus_id int(10) NOT NULL,
  analysis_id smallint(5) NOT NULL,
  score double(8,2) DEFAULT NULL
);

CREATE INDEX seq_region_idx14 ON repeat_feature (seq_region_id, seq_region_start);

CREATE INDEX repeat_idx ON repeat_feature (repeat_consensus_id);

CREATE INDEX analysis_idx10 ON repeat_feature (analysis_id);

--
-- Table: seq_region
--
CREATE TABLE seq_region (
  seq_region_id INTEGER PRIMARY KEY NOT NULL,
  name varchar(40) NOT NULL,
  coord_system_id int(10) NOT NULL,
  length int(10) NOT NULL
);

CREATE INDEX cs_idx ON seq_region (coord_system_id);

CREATE UNIQUE INDEX name_cs_idx ON seq_region (name, coord_system_id);

--
-- Table: seq_region_attrib
--
CREATE TABLE seq_region_attrib (
  seq_region_id int(10) NOT NULL DEFAULT 0,
  attrib_type_id smallint(5) NOT NULL DEFAULT 0,
  value text NOT NULL
);

CREATE INDEX type_val_idx03 ON seq_region_attrib (attrib_type_id, value);

CREATE INDEX val_only_idx03 ON seq_region_attrib (value);

CREATE INDEX seq_region_idx15 ON seq_region_attrib (seq_region_id);

--
-- Table: seq_region_mapping
--
CREATE TABLE seq_region_mapping (
  external_seq_region_id int(10) NOT NULL,
  internal_seq_region_id int(10) NOT NULL,
  mapping_set_id int(10) NOT NULL
);

CREATE INDEX mapping_set_idx ON seq_region_mapping (mapping_set_id);

--
-- Table: seq_region_synonym
--
CREATE TABLE seq_region_synonym (
  seq_region_synonym_id INTEGER PRIMARY KEY NOT NULL,
  seq_region_id int(10) NOT NULL,
  synonym varchar(50) NOT NULL,
  external_db_id int(10) DEFAULT NULL
);

CREATE UNIQUE INDEX syn_idx ON seq_region_synonym (synonym);

--
-- Table: simple_feature
--
CREATE TABLE simple_feature (
  simple_feature_id INTEGER PRIMARY KEY NOT NULL,
  seq_region_id int(10) NOT NULL,
  seq_region_start int(10) NOT NULL,
  seq_region_end int(10) NOT NULL,
  seq_region_strand tinyint(1) NOT NULL,
  display_label varchar(255) NOT NULL,
  analysis_id smallint(5) NOT NULL,
  score double(8,2) DEFAULT NULL
);

CREATE INDEX seq_region_idx16 ON simple_feature (seq_region_id, seq_region_start);

CREATE INDEX analysis_idx11 ON simple_feature (analysis_id);

CREATE INDEX hit_idx03 ON simple_feature (display_label);

--
-- Table: stable_id_event
--
CREATE TABLE stable_id_event (
  old_stable_id varchar(128) DEFAULT NULL,
  old_version smallint(6) DEFAULT NULL,
  new_stable_id varchar(128) DEFAULT NULL,
  new_version smallint(6) DEFAULT NULL,
  mapping_session_id int(10) NOT NULL DEFAULT 0,
  type enum(11) NOT NULL,
  score float(8,2) NOT NULL DEFAULT 0
);

CREATE INDEX new_idx ON stable_id_event (new_stable_id);

CREATE INDEX old_idx ON stable_id_event (old_stable_id);

CREATE UNIQUE INDEX uni_idx ON stable_id_event (mapping_session_id, old_stable_id, new_stable_id, type);

--
-- Table: supporting_feature
--
CREATE TABLE supporting_feature (
  exon_id int(10) NOT NULL DEFAULT 0,
  feature_type enum(21) DEFAULT NULL,
  feature_id int(10) NOT NULL DEFAULT 0
);

CREATE INDEX feature_idx ON supporting_feature (feature_type, feature_id);

CREATE UNIQUE INDEX all_idx02 ON supporting_feature (exon_id, feature_type, feature_id);

--
-- Table: transcript
--
CREATE TABLE transcript (
  transcript_id INTEGER PRIMARY KEY NOT NULL,
  gene_id int(10) DEFAULT NULL,
  analysis_id smallint(5) NOT NULL,
  seq_region_id int(10) NOT NULL,
  seq_region_start int(10) NOT NULL,
  seq_region_end int(10) NOT NULL,
  seq_region_strand tinyint(2) NOT NULL,
  display_xref_id int(10) DEFAULT NULL,
  source varchar(40) NOT NULL DEFAULT 'ensembl',
  biotype varchar(40) NOT NULL,
  status enum(19) DEFAULT NULL,
  description text,
  is_current tinyint(1) NOT NULL DEFAULT 1,
  canonical_translation_id int(10) DEFAULT NULL,
  stable_id varchar(128) DEFAULT NULL,
  version smallint(5) NOT NULL DEFAULT 1,
  created_date datetime NOT NULL DEFAULT '0000-00-00 00:00:00',
  modified_date datetime NOT NULL DEFAULT '0000-00-00 00:00:00'
);

CREATE INDEX seq_region_idx17 ON transcript (seq_region_id, seq_region_start);

CREATE INDEX gene_index ON transcript (gene_id);

CREATE INDEX xref_id_index02 ON transcript (display_xref_id);

CREATE INDEX analysis_idx12 ON transcript (analysis_id);

CREATE INDEX stable_id_idx05 ON transcript (stable_id, version);

CREATE UNIQUE INDEX canonical_translation_idx ON transcript (canonical_translation_id);

--
-- Table: transcript_attrib
--
CREATE TABLE transcript_attrib (
  transcript_id int(10) NOT NULL DEFAULT 0,
  attrib_type_id smallint(5) NOT NULL DEFAULT 0,
  value text NOT NULL
);

CREATE INDEX type_val_idx04 ON transcript_attrib (attrib_type_id, value);

CREATE INDEX val_only_idx04 ON transcript_attrib (value);

CREATE INDEX transcript_idx03 ON transcript_attrib (transcript_id);

--
-- Table: transcript_intron_supporting_evidence
--
CREATE TABLE transcript_intron_supporting_evidence (
  transcript_id int(10) NOT NULL,
  intron_supporting_evidence_id int(10) NOT NULL,
  previous_exon_id int(10) NOT NULL,
  next_exon_id int(10) NOT NULL,
  PRIMARY KEY (intron_supporting_evidence_id, transcript_id)
);

CREATE INDEX transcript_idx04 ON transcript_intron_supporting_evidence (transcript_id);

--
-- Table: transcript_supporting_feature
--
CREATE TABLE transcript_supporting_feature (
  transcript_id int(10) NOT NULL DEFAULT 0,
  feature_type enum(21) DEFAULT NULL,
  feature_id int(10) NOT NULL DEFAULT 0
);

CREATE INDEX feature_idx02 ON transcript_supporting_feature (feature_type, feature_id);

CREATE UNIQUE INDEX all_idx03 ON transcript_supporting_feature (transcript_id, feature_type, feature_id);

--
-- Table: translation
--
CREATE TABLE translation (
  translation_id INTEGER PRIMARY KEY NOT NULL,
  transcript_id int(10) NOT NULL,
  seq_start int(10) NOT NULL,
  start_exon_id int(10) NOT NULL,
  seq_end int(10) NOT NULL,
  end_exon_id int(10) NOT NULL,
  stable_id varchar(128) DEFAULT NULL,
  version smallint(5) NOT NULL DEFAULT 1,
  created_date datetime NOT NULL DEFAULT '0000-00-00 00:00:00',
  modified_date datetime NOT NULL DEFAULT '0000-00-00 00:00:00'
);

CREATE INDEX transcript_idx05 ON translation (transcript_id);

CREATE INDEX stable_id_idx06 ON translation (stable_id, version);

--
-- Table: translation_attrib
--
CREATE TABLE translation_attrib (
  translation_id int(10) NOT NULL DEFAULT 0,
  attrib_type_id smallint(5) NOT NULL DEFAULT 0,
  value text NOT NULL
);

CREATE INDEX type_val_idx05 ON translation_attrib (attrib_type_id, value);

CREATE INDEX val_only_idx05 ON translation_attrib (value);

CREATE INDEX translation_idx03 ON translation_attrib (translation_id);

--
-- Table: unmapped_object
--
CREATE TABLE unmapped_object (
  unmapped_object_id INTEGER PRIMARY KEY NOT NULL,
  type enum(6) NOT NULL,
  analysis_id smallint(5) NOT NULL,
  external_db_id int(10) DEFAULT NULL,
  identifier varchar(255) NOT NULL,
  unmapped_reason_id int(10) NOT NULL,
  query_score double(8,2) DEFAULT NULL,
  target_score double(8,2) DEFAULT NULL,
  ensembl_id int(10) DEFAULT 0,
  ensembl_object_type enum(11) DEFAULT 'RawContig',
  parent varchar(255) DEFAULT NULL
);

CREATE INDEX id_idx02 ON unmapped_object (identifier);

CREATE INDEX anal_exdb_idx ON unmapped_object (analysis_id, external_db_id);

--
-- Table: unmapped_reason
--
CREATE TABLE unmapped_reason (
  unmapped_reason_id INTEGER PRIMARY KEY NOT NULL,
  summary_description varchar(255) DEFAULT NULL,
  full_description varchar(255) DEFAULT NULL
);

--
-- Table: xref
--
CREATE TABLE xref (
  xref_id INTEGER PRIMARY KEY NOT NULL,
  external_db_id int(10) DEFAULT NULL,
  dbprimary_acc varchar(50) NOT NULL,
  display_label varchar(128) NOT NULL,
  version varchar(10) NOT NULL DEFAULT '0',
  description text,
  info_type enum(18) NOT NULL DEFAULT 'NONE',
  info_text varchar(255) NOT NULL DEFAULT ''
);

CREATE INDEX display_index ON xref (display_label);

CREATE INDEX info_type_idx ON xref (info_type);

CREATE UNIQUE INDEX id_index ON xref (dbprimary_acc, external_db_id, info_type, info_text, version);

COMMIT;
