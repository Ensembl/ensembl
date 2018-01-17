-- 
-- Created by SQL::Translator::Producer::SQLite
-- Created on Fri Jan 12 13:37:38 2018
-- 

BEGIN TRANSACTION;

--
-- Table: alt_allele
--
CREATE TABLE alt_allele (
  alt_allele_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  alt_allele_group_id integer NOT NULL,
  gene_id integer NOT NULL
);

CREATE UNIQUE INDEX gene_idx ON alt_allele (gene_id);

--
-- Table: alt_allele_attrib
--
CREATE TABLE alt_allele_attrib (
  alt_allele_id integer,
  attrib enum
);

--
-- Table: alt_allele_group
--
CREATE TABLE alt_allele_group (
  alt_allele_group_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL
);

--
-- Table: analysis
--
CREATE TABLE analysis (
  analysis_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  created datetime,
  logic_name varchar(40) NOT NULL DEFAULT '',
  db varchar(120),
  db_version varchar(40),
  db_file varchar(120),
  program varchar(80),
  program_version varchar(40),
  program_file varchar(80),
  parameters text,
  module varchar(80),
  module_version varchar(40),
  gff_source varchar(40),
  gff_feature varchar(40)
);

CREATE UNIQUE INDEX logic_name_idx ON analysis (logic_name);

--
-- Table: analysis_description
--
CREATE TABLE analysis_description (
  analysis_id integer NOT NULL DEFAULT 0,
  description text,
  display_label varchar(255),
  displayable tinyint NOT NULL DEFAULT 1,
  web_data text
);

CREATE UNIQUE INDEX analysis_idx ON analysis_description (analysis_id);

--
-- Table: assembly
--
CREATE TABLE assembly (
  asm_seq_region_id integer NOT NULL DEFAULT 0,
  cmp_seq_region_id integer NOT NULL DEFAULT 0,
  asm_start integer NOT NULL DEFAULT 0,
  asm_end integer NOT NULL DEFAULT 0,
  cmp_start integer NOT NULL DEFAULT 0,
  cmp_end integer NOT NULL DEFAULT 0,
  ori tinyint NOT NULL DEFAULT 0
);

CREATE UNIQUE INDEX all_idx ON assembly (asm_seq_region_id, cmp_seq_region_id, asm_start, asm_end, cmp_start, cmp_end, ori);

--
-- Table: assembly_exception
--
CREATE TABLE assembly_exception (
  assembly_exception_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  seq_region_id integer NOT NULL DEFAULT 0,
  seq_region_start integer NOT NULL DEFAULT 0,
  seq_region_end integer NOT NULL DEFAULT 0,
  exc_type enum NOT NULL DEFAULT 'HAP',
  exc_seq_region_id integer NOT NULL DEFAULT 0,
  exc_seq_region_start integer NOT NULL DEFAULT 0,
  exc_seq_region_end integer NOT NULL DEFAULT 0,
  ori integer NOT NULL DEFAULT 0
);

--
-- Table: associated_group
--
CREATE TABLE associated_group (
  associated_group_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  description varchar(128)
);

--
-- Table: associated_xref
--
CREATE TABLE associated_xref (
  associated_xref_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  object_xref_id integer NOT NULL DEFAULT 0,
  xref_id integer NOT NULL DEFAULT 0,
  source_xref_id integer,
  condition_type varchar(128),
  associated_group_id integer,
  rank integer DEFAULT 0
);

CREATE UNIQUE INDEX object_associated_source_type_idx ON associated_xref (object_xref_id, xref_id, source_xref_id, condition_type, associated_group_id);

--
-- Table: attrib_type
--
CREATE TABLE attrib_type (
  attrib_type_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  code varchar(20) NOT NULL DEFAULT '',
  name varchar(255) NOT NULL DEFAULT '',
  description text
);

CREATE UNIQUE INDEX code_idx ON attrib_type (code);

--
-- Table: coord_system
--
CREATE TABLE coord_system (
  coord_system_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  species_id integer NOT NULL DEFAULT 1,
  name varchar(40) NOT NULL,
  version varchar(255),
  rank integer NOT NULL,
  attrib varchar
);

CREATE UNIQUE INDEX name_idx ON coord_system (name, version, species_id);

CREATE UNIQUE INDEX rank_idx ON coord_system (rank, species_id);

--
-- Table: data_file
--
CREATE TABLE data_file (
  data_file_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  coord_system_id integer NOT NULL,
  analysis_id smallint NOT NULL,
  name varchar(100) NOT NULL,
  version_lock tinyint NOT NULL DEFAULT 0,
  absolute tinyint NOT NULL DEFAULT 0,
  url text,
  file_type enum
);

CREATE UNIQUE INDEX df_unq_idx ON data_file (coord_system_id, analysis_id, name, file_type);

--
-- Table: density_feature
--
CREATE TABLE density_feature (
  density_feature_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  density_type_id integer NOT NULL DEFAULT 0,
  seq_region_id integer NOT NULL DEFAULT 0,
  seq_region_start integer NOT NULL DEFAULT 0,
  seq_region_end integer NOT NULL DEFAULT 0,
  density_value float NOT NULL DEFAULT 0
);

--
-- Table: density_type
--
CREATE TABLE density_type (
  density_type_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  analysis_id integer NOT NULL DEFAULT 0,
  block_size integer NOT NULL DEFAULT 0,
  region_features integer NOT NULL DEFAULT 0,
  value_type enum NOT NULL DEFAULT 'sum'
);

CREATE UNIQUE INDEX analysis_id ON density_type (analysis_id, block_size, region_features);

--
-- Table: dependent_xref
--
CREATE TABLE dependent_xref (
  object_xref_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  master_xref_id integer NOT NULL,
  dependent_xref_id integer NOT NULL
);

--
-- Table: ditag
--
CREATE TABLE ditag (
  ditag_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  name varchar(30),
  type varchar(30),
  tag_count smallint DEFAULT 1,
  sequence text
);

--
-- Table: ditag_feature
--
CREATE TABLE ditag_feature (
  ditag_feature_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  ditag_id integer NOT NULL DEFAULT 0,
  ditag_pair_id integer NOT NULL DEFAULT 0,
  seq_region_id integer NOT NULL DEFAULT 0,
  seq_region_start integer NOT NULL DEFAULT 0,
  seq_region_end integer NOT NULL DEFAULT 0,
  seq_region_strand tinyint NOT NULL DEFAULT 0,
  analysis_id integer NOT NULL DEFAULT 0,
  hit_start integer NOT NULL DEFAULT 0,
  hit_end integer NOT NULL DEFAULT 0,
  hit_strand tinyint NOT NULL DEFAULT 0,
  cigar_line text,
  ditag_side char(1) DEFAULT ''
);

--
-- Table: dna
--
CREATE TABLE dna (
  seq_region_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL DEFAULT 0,
  sequence mediumtext NOT NULL
);

--
-- Table: dna_align_feature
--
CREATE TABLE dna_align_feature (
  dna_align_feature_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  seq_region_id integer NOT NULL DEFAULT 0,
  seq_region_start integer NOT NULL DEFAULT 0,
  seq_region_end integer NOT NULL DEFAULT 0,
  seq_region_strand tinyint NOT NULL DEFAULT 0,
  hit_start integer NOT NULL DEFAULT 0,
  hit_end integer NOT NULL DEFAULT 0,
  hit_strand tinyint NOT NULL DEFAULT 0,
  hit_name varchar(40) NOT NULL DEFAULT '',
  analysis_id integer NOT NULL DEFAULT 0,
  score double precision,
  evalue double precision,
  perc_ident float,
  cigar_line text,
  external_db_id smallint,
  hcoverage double precision,
  align_type enum DEFAULT 'ensembl'
);

--
-- Table: dna_align_feature_attrib
--
CREATE TABLE dna_align_feature_attrib (
  dna_align_feature_id integer NOT NULL,
  attrib_type_id smallint NOT NULL,
  value text NOT NULL
);

CREATE UNIQUE INDEX dna_align_feature_attribx ON dna_align_feature_attrib (dna_align_feature_id, attrib_type_id, value);

--
-- Table: exon
--
CREATE TABLE exon (
  exon_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  seq_region_id integer NOT NULL,
  seq_region_start integer NOT NULL,
  seq_region_end integer NOT NULL,
  seq_region_strand tinyint NOT NULL,
  phase tinyint NOT NULL,
  end_phase tinyint NOT NULL,
  is_current tinyint NOT NULL DEFAULT 1,
  is_constitutive tinyint NOT NULL DEFAULT 0,
  stable_id varchar(128),
  version smallint,
  created_date datetime,
  modified_date datetime
);

--
-- Table: exon_transcript
--
CREATE TABLE exon_transcript (
  exon_id integer NOT NULL DEFAULT 0,
  transcript_id integer NOT NULL DEFAULT 0,
  rank integer NOT NULL DEFAULT 0,
  PRIMARY KEY (exon_id, transcript_id, rank)
);

--
-- Table: external_db
--
CREATE TABLE external_db (
  external_db_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL DEFAULT 0,
  db_name varchar(27) NOT NULL DEFAULT '',
  db_release varchar(40) NOT NULL DEFAULT '',
  status enum NOT NULL DEFAULT 'KNOWNXREF',
  priority integer NOT NULL DEFAULT 0,
  db_display_name varchar(255),
  type enum,
  secondary_db_name varchar(255),
  secondary_db_table varchar(255),
  description text
);

CREATE UNIQUE INDEX db_name_db_release_idx ON external_db (db_name, db_release);

--
-- Table: external_synonym
--
CREATE TABLE external_synonym (
  xref_id integer NOT NULL DEFAULT 0,
  synonym varchar(40) NOT NULL DEFAULT '',
  PRIMARY KEY (xref_id, synonym)
);

--
-- Table: gene
--
CREATE TABLE gene (
  gene_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  biotype varchar(40) NOT NULL,
  analysis_id smallint NOT NULL,
  seq_region_id integer NOT NULL,
  seq_region_start integer NOT NULL,
  seq_region_end integer NOT NULL,
  seq_region_strand tinyint NOT NULL,
  display_xref_id integer,
  source varchar(40) NOT NULL,
  description text,
  is_current tinyint NOT NULL DEFAULT 1,
  canonical_transcript_id integer NOT NULL,
  stable_id varchar(128),
  version smallint,
  created_date datetime,
  modified_date datetime
);

--
-- Table: gene_archive
--
CREATE TABLE gene_archive (
  gene_stable_id varchar(128) NOT NULL DEFAULT '',
  gene_version smallint NOT NULL DEFAULT 0,
  transcript_stable_id varchar(128) NOT NULL DEFAULT '',
  transcript_version smallint NOT NULL DEFAULT 0,
  translation_stable_id varchar(128) NOT NULL DEFAULT '',
  translation_version smallint NOT NULL DEFAULT 0,
  peptide_archive_id integer NOT NULL DEFAULT 0,
  mapping_session_id integer NOT NULL DEFAULT 0
);

--
-- Table: gene_attrib
--
CREATE TABLE gene_attrib (
  gene_id integer NOT NULL DEFAULT 0,
  attrib_type_id smallint NOT NULL DEFAULT 0,
  value text NOT NULL
);

CREATE UNIQUE INDEX gene_attribx ON gene_attrib (gene_id, attrib_type_id, value);

--
-- Table: genome_statistics
--
CREATE TABLE genome_statistics (
  genome_statistics_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  statistic varchar(128) NOT NULL,
  value bigint NOT NULL DEFAULT 0,
  species_id integer DEFAULT 1,
  attrib_type_id integer,
  timestamp datetime
);

CREATE UNIQUE INDEX stats_uniq ON genome_statistics (statistic, attrib_type_id, species_id);

--
-- Table: identity_xref
--
CREATE TABLE identity_xref (
  object_xref_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL DEFAULT 0,
  xref_identity integer,
  ensembl_identity integer,
  xref_start integer,
  xref_end integer,
  ensembl_start integer,
  ensembl_end integer,
  cigar_line text,
  score double precision,
  evalue double precision
);

--
-- Table: interpro
--
CREATE TABLE interpro (
  interpro_ac varchar(40) NOT NULL DEFAULT '',
  id varchar(40) NOT NULL DEFAULT ''
);

CREATE UNIQUE INDEX accession_idx ON interpro (interpro_ac, id);

--
-- Table: intron_supporting_evidence
--
CREATE TABLE intron_supporting_evidence (
  intron_supporting_evidence_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  analysis_id smallint NOT NULL,
  seq_region_id integer NOT NULL,
  seq_region_start integer NOT NULL,
  seq_region_end integer NOT NULL,
  seq_region_strand tinyint NOT NULL,
  hit_name varchar(100) NOT NULL,
  score decimal(10,3),
  score_type enum DEFAULT 'NONE',
  is_splice_canonical tinyint NOT NULL DEFAULT 0
);

CREATE UNIQUE INDEX analysis_id02 ON intron_supporting_evidence (analysis_id, seq_region_id, seq_region_start, seq_region_end, seq_region_strand, hit_name);

--
-- Table: karyotype
--
CREATE TABLE karyotype (
  karyotype_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  seq_region_id integer NOT NULL DEFAULT 0,
  seq_region_start integer NOT NULL DEFAULT 0,
  seq_region_end integer NOT NULL DEFAULT 0,
  band varchar(40),
  stain varchar(40)
);

--
-- Table: map
--
CREATE TABLE map (
  map_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  map_name varchar(30) NOT NULL DEFAULT ''
);

--
-- Table: mapping_session
--
CREATE TABLE mapping_session (
  mapping_session_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  old_db_name varchar(80) NOT NULL DEFAULT '',
  new_db_name varchar(80) NOT NULL DEFAULT '',
  old_release varchar(5) NOT NULL DEFAULT '',
  new_release varchar(5) NOT NULL DEFAULT '',
  old_assembly varchar(20) NOT NULL DEFAULT '',
  new_assembly varchar(20) NOT NULL DEFAULT '',
  created datetime
);

--
-- Table: mapping_set
--
CREATE TABLE mapping_set (
  mapping_set_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  internal_schema_build varchar(20) NOT NULL,
  external_schema_build varchar(20) NOT NULL
);

CREATE UNIQUE INDEX mapping_idx ON mapping_set (internal_schema_build, external_schema_build);

--
-- Table: marker
--
CREATE TABLE marker (
  marker_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  display_marker_synonym_id integer,
  left_primer varchar(100) NOT NULL DEFAULT '',
  right_primer varchar(100) NOT NULL DEFAULT '',
  min_primer_dist integer NOT NULL DEFAULT 0,
  max_primer_dist integer NOT NULL DEFAULT 0,
  priority integer,
  type enum
);

--
-- Table: marker_feature
--
CREATE TABLE marker_feature (
  marker_feature_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  marker_id integer NOT NULL DEFAULT 0,
  seq_region_id integer NOT NULL DEFAULT 0,
  seq_region_start integer NOT NULL DEFAULT 0,
  seq_region_end integer NOT NULL DEFAULT 0,
  analysis_id integer NOT NULL DEFAULT 0,
  map_weight integer
);

--
-- Table: marker_map_location
--
CREATE TABLE marker_map_location (
  marker_id integer NOT NULL DEFAULT 0,
  map_id integer NOT NULL DEFAULT 0,
  chromosome_name varchar(15) NOT NULL DEFAULT '',
  marker_synonym_id integer NOT NULL DEFAULT 0,
  position varchar(15) NOT NULL DEFAULT '',
  lod_score double precision,
  PRIMARY KEY (marker_id, map_id)
);

--
-- Table: marker_synonym
--
CREATE TABLE marker_synonym (
  marker_synonym_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  marker_id integer NOT NULL DEFAULT 0,
  source varchar(20),
  name varchar(30)
);

--
-- Table: meta
--
CREATE TABLE meta (
  meta_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  species_id integer DEFAULT 1,
  meta_key varchar(40) NOT NULL,
  meta_value varchar(255) NOT NULL
);

CREATE UNIQUE INDEX species_key_value_idx ON meta (species_id, meta_key, meta_value);

--
-- Table: meta_coord
--
CREATE TABLE meta_coord (
  table_name varchar(40) NOT NULL DEFAULT '',
  coord_system_id integer NOT NULL DEFAULT 0,
  max_length integer
);

CREATE UNIQUE INDEX cs_table_name_idx ON meta_coord (coord_system_id, table_name);

--
-- Table: misc_attrib
--
CREATE TABLE misc_attrib (
  misc_feature_id integer NOT NULL DEFAULT 0,
  attrib_type_id smallint NOT NULL DEFAULT 0,
  value text NOT NULL
);

CREATE UNIQUE INDEX misc_attribx ON misc_attrib (misc_feature_id, attrib_type_id, value);

--
-- Table: misc_feature
--
CREATE TABLE misc_feature (
  misc_feature_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  seq_region_id integer NOT NULL DEFAULT 0,
  seq_region_start integer NOT NULL DEFAULT 0,
  seq_region_end integer NOT NULL DEFAULT 0,
  seq_region_strand tinyint NOT NULL DEFAULT 0
);

--
-- Table: misc_feature_misc_set
--
CREATE TABLE misc_feature_misc_set (
  misc_feature_id integer NOT NULL DEFAULT 0,
  misc_set_id smallint NOT NULL DEFAULT 0,
  PRIMARY KEY (misc_feature_id, misc_set_id)
);

--
-- Table: misc_set
--
CREATE TABLE misc_set (
  misc_set_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  code varchar(25) NOT NULL DEFAULT '',
  name varchar(255) NOT NULL DEFAULT '',
  description text NOT NULL,
  max_length integer NOT NULL DEFAULT 0
);

CREATE UNIQUE INDEX code_idx02 ON misc_set (code);

--
-- Table: object_xref
--
CREATE TABLE object_xref (
  object_xref_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  ensembl_id integer NOT NULL DEFAULT 0,
  ensembl_object_type enum NOT NULL DEFAULT 'RawContig',
  xref_id integer NOT NULL,
  linkage_annotation varchar(255),
  analysis_id smallint NOT NULL
);

CREATE UNIQUE INDEX xref_idx ON object_xref (xref_id, ensembl_object_type, ensembl_id, analysis_id);

--
-- Table: ontology_xref
--
CREATE TABLE ontology_xref (
  object_xref_id integer NOT NULL DEFAULT 0,
  linkage_type varchar(3),
  source_xref_id integer
);

CREATE UNIQUE INDEX object_source_type_idx ON ontology_xref (object_xref_id, source_xref_id, linkage_type);

--
-- Table: operon
--
CREATE TABLE operon (
  operon_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  seq_region_id integer NOT NULL,
  seq_region_start integer NOT NULL,
  seq_region_end integer NOT NULL,
  seq_region_strand tinyint NOT NULL,
  display_label varchar(255),
  analysis_id smallint NOT NULL,
  stable_id varchar(128),
  version smallint,
  created_date datetime,
  modified_date datetime
);

--
-- Table: operon_transcript
--
CREATE TABLE operon_transcript (
  operon_transcript_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  seq_region_id integer NOT NULL,
  seq_region_start integer NOT NULL,
  seq_region_end integer NOT NULL,
  seq_region_strand tinyint NOT NULL,
  operon_id integer NOT NULL,
  display_label varchar(255),
  analysis_id smallint NOT NULL,
  stable_id varchar(128),
  version smallint,
  created_date datetime,
  modified_date datetime
);

--
-- Table: operon_transcript_gene
--
CREATE TABLE operon_transcript_gene (
  operon_transcript_id integer,
  gene_id integer
);

--
-- Table: peptide_archive
--
CREATE TABLE peptide_archive (
  peptide_archive_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  md5_checksum varchar(32),
  peptide_seq mediumtext NOT NULL
);

--
-- Table: prediction_exon
--
CREATE TABLE prediction_exon (
  prediction_exon_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  prediction_transcript_id integer NOT NULL DEFAULT 0,
  exon_rank smallint NOT NULL DEFAULT 0,
  seq_region_id integer NOT NULL DEFAULT 0,
  seq_region_start integer NOT NULL DEFAULT 0,
  seq_region_end integer NOT NULL DEFAULT 0,
  seq_region_strand tinyint NOT NULL DEFAULT 0,
  start_phase tinyint NOT NULL DEFAULT 0,
  score double precision,
  p_value double precision
);

--
-- Table: prediction_transcript
--
CREATE TABLE prediction_transcript (
  prediction_transcript_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  seq_region_id integer NOT NULL DEFAULT 0,
  seq_region_start integer NOT NULL DEFAULT 0,
  seq_region_end integer NOT NULL DEFAULT 0,
  seq_region_strand tinyint NOT NULL DEFAULT 0,
  analysis_id integer,
  display_label varchar(255)
);

--
-- Table: protein_align_feature
--
CREATE TABLE protein_align_feature (
  protein_align_feature_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  seq_region_id integer NOT NULL DEFAULT 0,
  seq_region_start integer NOT NULL DEFAULT 0,
  seq_region_end integer NOT NULL DEFAULT 0,
  seq_region_strand tinyint NOT NULL DEFAULT 1,
  hit_start integer NOT NULL DEFAULT 0,
  hit_end integer NOT NULL DEFAULT 0,
  hit_name varchar(40) NOT NULL DEFAULT '',
  analysis_id integer NOT NULL DEFAULT 0,
  score double precision,
  evalue double precision,
  perc_ident float,
  cigar_line text,
  external_db_id smallint,
  hcoverage double precision,
  align_type enum DEFAULT 'ensembl'
);

--
-- Table: protein_feature
--
CREATE TABLE protein_feature (
  protein_feature_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  translation_id integer NOT NULL DEFAULT 0,
  seq_start integer NOT NULL DEFAULT 0,
  seq_end integer NOT NULL DEFAULT 0,
  hit_start integer NOT NULL DEFAULT 0,
  hit_end integer NOT NULL DEFAULT 0,
  hit_name varchar(40) NOT NULL,
  analysis_id integer NOT NULL DEFAULT 0,
  score double precision NOT NULL DEFAULT 0,
  evalue double precision,
  perc_ident float,
  external_data text,
  hit_description text,
  cigar_line text,
  align_type enum
);

CREATE UNIQUE INDEX aln_idx ON protein_feature (translation_id, hit_name, seq_start, seq_end, hit_start, hit_end, analysis_id);

--
-- Table: repeat_consensus
--
CREATE TABLE repeat_consensus (
  repeat_consensus_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  repeat_name varchar(255) NOT NULL DEFAULT '',
  repeat_class varchar(100) NOT NULL DEFAULT '',
  repeat_type varchar(40) NOT NULL DEFAULT '',
  repeat_consensus text
);

--
-- Table: repeat_feature
--
CREATE TABLE repeat_feature (
  repeat_feature_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  seq_region_id integer NOT NULL DEFAULT 0,
  seq_region_start integer NOT NULL DEFAULT 0,
  seq_region_end integer NOT NULL DEFAULT 0,
  seq_region_strand tinyint NOT NULL DEFAULT 1,
  repeat_start integer NOT NULL DEFAULT 0,
  repeat_end integer NOT NULL DEFAULT 0,
  repeat_consensus_id integer NOT NULL DEFAULT 0,
  analysis_id integer NOT NULL DEFAULT 0,
  score double precision
);

--
-- Table: seq_region
--
CREATE TABLE seq_region (
  seq_region_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  name varchar(255) NOT NULL,
  coord_system_id integer NOT NULL DEFAULT 0,
  length integer NOT NULL DEFAULT 0
);

CREATE UNIQUE INDEX name_cs_idx ON seq_region (name, coord_system_id);

--
-- Table: seq_region_attrib
--
CREATE TABLE seq_region_attrib (
  seq_region_id integer NOT NULL DEFAULT 0,
  attrib_type_id smallint NOT NULL DEFAULT 0,
  value text NOT NULL
);

CREATE UNIQUE INDEX region_attribx ON seq_region_attrib (seq_region_id, attrib_type_id, value);

--
-- Table: seq_region_mapping
--
CREATE TABLE seq_region_mapping (
  external_seq_region_id integer NOT NULL,
  internal_seq_region_id integer NOT NULL,
  mapping_set_id integer NOT NULL
);

--
-- Table: seq_region_synonym
--
CREATE TABLE seq_region_synonym (
  seq_region_synonym_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  seq_region_id integer NOT NULL,
  synonym varchar(250) NOT NULL,
  external_db_id smallint
);

CREATE UNIQUE INDEX syn_idx ON seq_region_synonym (synonym, seq_region_id);

--
-- Table: simple_feature
--
CREATE TABLE simple_feature (
  simple_feature_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  seq_region_id integer NOT NULL DEFAULT 0,
  seq_region_start integer NOT NULL DEFAULT 0,
  seq_region_end integer NOT NULL DEFAULT 0,
  seq_region_strand tinyint NOT NULL DEFAULT 0,
  display_label varchar(40) NOT NULL DEFAULT '',
  analysis_id integer NOT NULL DEFAULT 0,
  score double precision
);

--
-- Table: stable_id_event
--
CREATE TABLE stable_id_event (
  old_stable_id varchar(128),
  old_version smallint,
  new_stable_id varchar(128),
  new_version smallint,
  mapping_session_id integer NOT NULL DEFAULT 0,
  type enum NOT NULL DEFAULT 'gene',
  score float NOT NULL DEFAULT 0
);

CREATE UNIQUE INDEX uni_idx ON stable_id_event (mapping_session_id, old_stable_id, old_version, new_stable_id, new_version, type);

--
-- Table: supporting_feature
--
CREATE TABLE supporting_feature (
  exon_id integer NOT NULL DEFAULT 0,
  feature_type enum,
  feature_id integer NOT NULL DEFAULT 0
);

CREATE UNIQUE INDEX all_idx02 ON supporting_feature (exon_id, feature_type, feature_id);

--
-- Table: transcript
--
CREATE TABLE transcript (
  transcript_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  gene_id integer,
  analysis_id smallint NOT NULL,
  seq_region_id integer NOT NULL,
  seq_region_start integer NOT NULL,
  seq_region_end integer NOT NULL,
  seq_region_strand tinyint NOT NULL,
  display_xref_id integer,
  source varchar(40) NOT NULL DEFAULT 'ensembl',
  biotype varchar(40) NOT NULL,
  description text,
  is_current tinyint NOT NULL DEFAULT 1,
  canonical_translation_id integer,
  stable_id varchar(128),
  version smallint,
  created_date datetime,
  modified_date datetime
);

CREATE UNIQUE INDEX canonical_translation_idx ON transcript (canonical_translation_id);

--
-- Table: transcript_attrib
--
CREATE TABLE transcript_attrib (
  transcript_id integer NOT NULL DEFAULT 0,
  attrib_type_id smallint NOT NULL DEFAULT 0,
  value text NOT NULL
);

CREATE UNIQUE INDEX transcript_attribx ON transcript_attrib (transcript_id, attrib_type_id, value);

--
-- Table: transcript_intron_supporting_evidence
--
CREATE TABLE transcript_intron_supporting_evidence (
  transcript_id integer NOT NULL,
  intron_supporting_evidence_id integer NOT NULL,
  previous_exon_id integer NOT NULL,
  next_exon_id integer NOT NULL,
  PRIMARY KEY (intron_supporting_evidence_id, transcript_id)
);

--
-- Table: transcript_supporting_feature
--
CREATE TABLE transcript_supporting_feature (
  transcript_id integer NOT NULL DEFAULT 0,
  feature_type enum,
  feature_id integer NOT NULL DEFAULT 0
);

CREATE UNIQUE INDEX all_idx03 ON transcript_supporting_feature (transcript_id, feature_type, feature_id);

--
-- Table: translation
--
CREATE TABLE translation (
  translation_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  transcript_id integer NOT NULL,
  seq_start integer NOT NULL,
  start_exon_id integer NOT NULL,
  seq_end integer NOT NULL,
  end_exon_id integer NOT NULL,
  stable_id varchar(128),
  version smallint,
  created_date datetime,
  modified_date datetime
);

--
-- Table: translation_attrib
--
CREATE TABLE translation_attrib (
  translation_id integer NOT NULL DEFAULT 0,
  attrib_type_id smallint NOT NULL DEFAULT 0,
  value text NOT NULL
);

CREATE UNIQUE INDEX translation_attribx ON translation_attrib (translation_id, attrib_type_id, value);

--
-- Table: unmapped_object
--
CREATE TABLE unmapped_object (
  unmapped_object_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  type enum NOT NULL,
  analysis_id integer NOT NULL,
  external_db_id integer,
  identifier varchar(255) NOT NULL,
  unmapped_reason_id integer NOT NULL,
  query_score double precision,
  target_score double precision,
  ensembl_id integer DEFAULT 0,
  ensembl_object_type enum DEFAULT 'RawContig',
  parent varchar(255)
);

CREATE UNIQUE INDEX unique_unmapped_obj_idx ON unmapped_object (ensembl_id, ensembl_object_type, identifier, unmapped_reason_id, parent, external_db_id);

--
-- Table: unmapped_reason
--
CREATE TABLE unmapped_reason (
  unmapped_reason_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  summary_description varchar(255),
  full_description varchar(255)
);

--
-- Table: xref
--
CREATE TABLE xref (
  xref_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  external_db_id integer NOT NULL,
  dbprimary_acc varchar(512) NOT NULL,
  display_label varchar(512) NOT NULL,
  version varchar(10),
  description text,
  info_type enum NOT NULL DEFAULT 'NONE',
  info_text varchar(255) NOT NULL DEFAULT ''
);

CREATE UNIQUE INDEX id_index ON xref (dbprimary_acc, external_db_id, info_type, info_text, version);

COMMIT;
