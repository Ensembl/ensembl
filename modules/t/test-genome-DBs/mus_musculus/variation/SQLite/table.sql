-- 
-- Created by SQL::Translator::Producer::SQLite
-- Created on Thu Nov 12 12:46:48 2015
-- 

BEGIN TRANSACTION;

--
-- Table: allele
--
CREATE TABLE allele (
  allele_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  variation_id integer NOT NULL,
  subsnp_id integer,
  allele_code_id integer NOT NULL,
  population_id integer,
  frequency float,
  count integer,
  frequency_submitter_handle integer
);

--
-- Table: allele_code
--
CREATE TABLE allele_code (
  allele_code_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  allele varchar(100)
);

CREATE UNIQUE INDEX allele_idx ON allele_code (allele);

--
-- Table: associate_study
--
CREATE TABLE associate_study (
  study1_id integer NOT NULL,
  study2_id integer NOT NULL,
  PRIMARY KEY (study1_id, study2_id)
);

--
-- Table: attrib
--
CREATE TABLE attrib (
  attrib_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL DEFAULT 0,
  attrib_type_id smallint NOT NULL DEFAULT 0,
  value text NOT NULL
);

CREATE UNIQUE INDEX type_val_idx ON attrib (attrib_type_id, value);

--
-- Table: attrib_set
--
CREATE TABLE attrib_set (
  attrib_set_id integer NOT NULL DEFAULT 0,
  attrib_id integer NOT NULL DEFAULT 0
);

CREATE UNIQUE INDEX set_idx ON attrib_set (attrib_set_id, attrib_id);

--
-- Table: attrib_type
--
CREATE TABLE attrib_type (
  attrib_type_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL DEFAULT 0,
  code varchar(20) NOT NULL DEFAULT '',
  name varchar(255) NOT NULL DEFAULT '',
  description text
);

CREATE UNIQUE INDEX code_idx ON attrib_type (code);

--
-- Table: compressed_genotype_region
--
CREATE TABLE compressed_genotype_region (
  sample_id integer NOT NULL,
  seq_region_id integer NOT NULL,
  seq_region_start integer NOT NULL,
  seq_region_end integer NOT NULL,
  seq_region_strand tinyint NOT NULL,
  genotypes blob
);

--
-- Table: compressed_genotype_var
--
CREATE TABLE compressed_genotype_var (
  variation_id integer NOT NULL,
  subsnp_id integer,
  genotypes blob
);

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
-- Table: display_group
--
CREATE TABLE display_group (
  display_group_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  display_priority integer NOT NULL,
  display_name varchar(255) NOT NULL
);

CREATE UNIQUE INDEX display_name ON display_group (display_name);

CREATE UNIQUE INDEX display_priority ON display_group (display_priority);

--
-- Table: failed_allele
--
CREATE TABLE failed_allele (
  failed_allele_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  allele_id integer NOT NULL,
  failed_description_id integer NOT NULL
);

CREATE UNIQUE INDEX allele_idx02 ON failed_allele (allele_id, failed_description_id);

--
-- Table: failed_description
--
CREATE TABLE failed_description (
  failed_description_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  description text NOT NULL
);

--
-- Table: failed_structural_variation
--
CREATE TABLE failed_structural_variation (
  failed_structural_variation_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  structural_variation_id integer NOT NULL,
  failed_description_id integer NOT NULL
);

CREATE UNIQUE INDEX structural_variation_idx ON failed_structural_variation (structural_variation_id, failed_description_id);

--
-- Table: failed_variation
--
CREATE TABLE failed_variation (
  failed_variation_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  variation_id integer NOT NULL,
  failed_description_id integer NOT NULL
);

CREATE UNIQUE INDEX variation_idx ON failed_variation (variation_id, failed_description_id);

--
-- Table: genotype_code
--
CREATE TABLE genotype_code (
  genotype_code_id integer NOT NULL,
  allele_code_id integer NOT NULL,
  haplotype_id tinyint NOT NULL,
  phased tinyint
);

--
-- Table: individual
--
CREATE TABLE individual (
  individual_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  name varchar(255),
  description text,
  gender enum NOT NULL DEFAULT 'Unknown',
  father_individual_id integer,
  mother_individual_id integer,
  individual_type_id integer NOT NULL DEFAULT 0
);

--
-- Table: individual_synonym
--
CREATE TABLE individual_synonym (
  synonym_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  individual_id integer NOT NULL,
  source_id integer NOT NULL,
  name varchar(255)
);

--
-- Table: individual_type
--
CREATE TABLE individual_type (
  individual_type_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  name varchar(255) NOT NULL,
  description text
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
  table_name varchar(40) NOT NULL,
  coord_system_id integer NOT NULL,
  max_length integer
);

CREATE UNIQUE INDEX table_name ON meta_coord (table_name, coord_system_id);

--
-- Table: motif_feature_variation
--
CREATE TABLE motif_feature_variation (
  motif_feature_variation_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  variation_feature_id integer NOT NULL,
  feature_stable_id varchar(128),
  motif_feature_id integer NOT NULL,
  allele_string text,
  somatic tinyint NOT NULL DEFAULT 0,
  consequence_types varchar,
  motif_name varchar(60),
  motif_start integer,
  motif_end integer,
  motif_score_delta float,
  in_informative_position tinyint NOT NULL DEFAULT 0
);

--
-- Table: phenotype
--
CREATE TABLE phenotype (
  phenotype_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  stable_id varchar(255),
  name varchar(50),
  description varchar(255)
);

CREATE UNIQUE INDEX desc_idx ON phenotype (description);

--
-- Table: phenotype_feature
--
CREATE TABLE phenotype_feature (
  phenotype_feature_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  phenotype_id integer,
  source_id integer,
  study_id integer,
  type enum,
  object_id varchar(255),
  is_significant tinyint DEFAULT 1,
  seq_region_id integer,
  seq_region_start integer,
  seq_region_end integer,
  seq_region_strand tinyint
);

--
-- Table: phenotype_feature_attrib
--
CREATE TABLE phenotype_feature_attrib (
  phenotype_feature_id integer NOT NULL,
  attrib_type_id integer,
  value varchar(255)
);

--
-- Table: population
--
CREATE TABLE population (
  population_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  name varchar(255),
  size integer,
  description text,
  collection tinyint DEFAULT 0,
  freqs_from_gts tinyint,
  display enum DEFAULT 'UNDISPLAYABLE',
  display_group_id tinyint
);

--
-- Table: population_genotype
--
CREATE TABLE population_genotype (
  population_genotype_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  variation_id integer NOT NULL,
  subsnp_id integer,
  genotype_code_id integer,
  frequency float,
  population_id integer,
  count integer
);

--
-- Table: population_structure
--
CREATE TABLE population_structure (
  super_population_id integer NOT NULL,
  sub_population_id integer NOT NULL
);

CREATE UNIQUE INDEX super_population_idx ON population_structure (super_population_id, sub_population_id);

--
-- Table: population_synonym
--
CREATE TABLE population_synonym (
  synonym_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  population_id integer NOT NULL,
  source_id integer NOT NULL,
  name varchar(255)
);

--
-- Table: protein_function_predictions
--
CREATE TABLE protein_function_predictions (
  translation_md5_id integer NOT NULL,
  analysis_attrib_id integer NOT NULL,
  prediction_matrix mediumblob,
  PRIMARY KEY (translation_md5_id, analysis_attrib_id)
);

--
-- Table: protein_function_predictions_attrib
--
CREATE TABLE protein_function_predictions_attrib (
  translation_md5_id integer NOT NULL,
  analysis_attrib_id integer NOT NULL,
  attrib_type_id integer NOT NULL,
  position_values blob,
  PRIMARY KEY (translation_md5_id, analysis_attrib_id, attrib_type_id)
);

--
-- Table: publication
--
CREATE TABLE publication (
  publication_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  title varchar(255),
  authors varchar(255),
  pmid integer,
  pmcid varchar(255),
  year integer,
  doi varchar(50),
  ucsc_id varchar(50)
);

--
-- Table: read_coverage
--
CREATE TABLE read_coverage (
  seq_region_id integer NOT NULL,
  seq_region_start integer NOT NULL,
  seq_region_end integer NOT NULL,
  level tinyint NOT NULL,
  sample_id integer NOT NULL
);

--
-- Table: regulatory_feature_variation
--
CREATE TABLE regulatory_feature_variation (
  regulatory_feature_variation_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  variation_feature_id integer NOT NULL,
  feature_stable_id varchar(128),
  feature_type text,
  allele_string text,
  somatic tinyint NOT NULL DEFAULT 0,
  consequence_types varchar
);

--
-- Table: sample
--
CREATE TABLE sample (
  sample_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  individual_id integer NOT NULL,
  name varchar(255),
  description text,
  study_id integer,
  display enum DEFAULT 'UNDISPLAYABLE',
  has_coverage tinyint NOT NULL DEFAULT 0,
  variation_set_id varchar
);

--
-- Table: sample_genotype_multiple_bp
--
CREATE TABLE sample_genotype_multiple_bp (
  variation_id integer NOT NULL,
  subsnp_id integer,
  allele_1 varchar(100),
  allele_2 varchar(100),
  sample_id integer
);

--
-- Table: sample_population
--
CREATE TABLE sample_population (
  sample_id integer NOT NULL,
  population_id integer NOT NULL
);

--
-- Table: seq_region
--
CREATE TABLE seq_region (
  seq_region_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  name varchar(40) NOT NULL,
  coord_system_id integer NOT NULL
);

CREATE UNIQUE INDEX name_cs_idx ON seq_region (name, coord_system_id);

--
-- Table: source
--
CREATE TABLE source (
  source_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  name varchar(24) NOT NULL,
  version integer,
  description varchar(255),
  url varchar(255),
  type enum,
  somatic_status enum DEFAULT 'germline',
  data_types varchar
);

--
-- Table: structural_variation
--
CREATE TABLE structural_variation (
  structural_variation_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  variation_name varchar(255),
  alias varchar(255),
  source_id integer NOT NULL,
  study_id integer,
  class_attrib_id integer NOT NULL DEFAULT 0,
  clinical_significance varchar,
  validation_status enum,
  is_evidence tinyint DEFAULT 0,
  somatic tinyint NOT NULL DEFAULT 0,
  copy_number tinyint
);

CREATE UNIQUE INDEX variation_name ON structural_variation (variation_name);

--
-- Table: structural_variation_association
--
CREATE TABLE structural_variation_association (
  structural_variation_id integer NOT NULL,
  supporting_structural_variation_id integer NOT NULL,
  PRIMARY KEY (structural_variation_id, supporting_structural_variation_id)
);

--
-- Table: structural_variation_feature
--
CREATE TABLE structural_variation_feature (
  structural_variation_feature_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  seq_region_id integer NOT NULL,
  outer_start integer,
  seq_region_start integer NOT NULL,
  inner_start integer,
  inner_end integer,
  seq_region_end integer NOT NULL,
  outer_end integer,
  seq_region_strand tinyint NOT NULL,
  structural_variation_id integer NOT NULL,
  variation_name varchar(255),
  source_id integer NOT NULL,
  study_id integer,
  class_attrib_id integer NOT NULL DEFAULT 0,
  allele_string longtext,
  is_evidence tinyint NOT NULL DEFAULT 0,
  variation_set_id varchar NOT NULL DEFAULT '',
  somatic tinyint NOT NULL DEFAULT 0,
  breakpoint_order tinyint,
  length integer
);

--
-- Table: structural_variation_sample
--
CREATE TABLE structural_variation_sample (
  structural_variation_sample_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  structural_variation_id integer NOT NULL,
  sample_id integer
);

--
-- Table: study
--
CREATE TABLE study (
  study_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  source_id integer NOT NULL,
  name varchar(255),
  description text,
  url varchar(255),
  external_reference varchar(255),
  study_type varchar(255)
);

--
-- Table: submitter_handle
--
CREATE TABLE submitter_handle (
  handle_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  handle varchar(25)
);

CREATE UNIQUE INDEX handle ON submitter_handle (handle);

--
-- Table: subsnp_handle
--
CREATE TABLE subsnp_handle (
  subsnp_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  handle varchar(20)
);

--
-- Table: subsnp_map
--
CREATE TABLE subsnp_map (
  variation_id integer NOT NULL,
  subsnp_id integer
);

--
-- Table: tagged_variation_feature
--
CREATE TABLE tagged_variation_feature (
  variation_feature_id integer NOT NULL,
  tagged_variation_feature_id integer,
  population_id integer NOT NULL
);

--
-- Table: tmp_sample_genotype_single_bp
--
CREATE TABLE tmp_sample_genotype_single_bp (
  variation_id integer NOT NULL,
  subsnp_id integer,
  allele_1 char(1),
  allele_2 char(1),
  sample_id integer NOT NULL
);

--
-- Table: transcript_variation
--
CREATE TABLE transcript_variation (
  transcript_variation_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  variation_feature_id integer NOT NULL,
  feature_stable_id varchar(128),
  allele_string text,
  somatic tinyint NOT NULL DEFAULT 0,
  consequence_types varchar,
  cds_start integer,
  cds_end integer,
  cdna_start integer,
  cdna_end integer,
  translation_start integer,
  translation_end integer,
  distance_to_transcript integer,
  codon_allele_string text,
  pep_allele_string text,
  hgvs_genomic text,
  hgvs_transcript text,
  hgvs_protein text,
  polyphen_prediction enum,
  sift_prediction enum,
  polyphen_score float,
  sift_score float,
  display integer DEFAULT 1
);

--
-- Table: translation_md5
--
CREATE TABLE translation_md5 (
  translation_md5_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  translation_md5 char(32) NOT NULL
);

CREATE UNIQUE INDEX md5_idx ON translation_md5 (translation_md5);

--
-- Table: variation
--
CREATE TABLE variation (
  variation_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  source_id integer NOT NULL,
  name varchar(255),
  ancestral_allele varchar(255),
  flipped tinyint,
  class_attrib_id integer DEFAULT 0,
  somatic tinyint NOT NULL DEFAULT 0,
  minor_allele varchar(50),
  minor_allele_freq float,
  minor_allele_count integer,
  clinical_significance varchar,
  evidence_attribs varchar,
  display integer DEFAULT 1
);

CREATE UNIQUE INDEX name ON variation (name);

--
-- Table: variation_attrib
--
CREATE TABLE variation_attrib (
  variation_id integer NOT NULL,
  attrib_id integer,
  value varchar(255)
);

--
-- Table: variation_citation
--
CREATE TABLE variation_citation (
  variation_id integer NOT NULL,
  publication_id integer NOT NULL,
  PRIMARY KEY (variation_id, publication_id)
);

--
-- Table: variation_feature
--
CREATE TABLE variation_feature (
  variation_feature_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  seq_region_id integer NOT NULL,
  seq_region_start integer NOT NULL,
  seq_region_end integer NOT NULL,
  seq_region_strand tinyint NOT NULL,
  variation_id integer NOT NULL,
  allele_string varchar(100),
  variation_name varchar(255),
  map_weight integer NOT NULL,
  flags varchar,
  source_id integer NOT NULL,
  consequence_types varchar NOT NULL DEFAULT 'intergenic_variant',
  variation_set_id varchar NOT NULL DEFAULT '',
  class_attrib_id integer DEFAULT 0,
  somatic tinyint NOT NULL DEFAULT 0,
  minor_allele varchar(50),
  minor_allele_freq float,
  minor_allele_count integer,
  alignment_quality double precision,
  evidence_attribs varchar,
  clinical_significance varchar,
  display integer DEFAULT 1
);

--
-- Table: variation_genename
--
CREATE TABLE variation_genename (
  variation_id integer NOT NULL,
  gene_name varchar(255) NOT NULL,
  PRIMARY KEY (variation_id, gene_name)
);

--
-- Table: variation_hgvs
--
CREATE TABLE variation_hgvs (
  variation_id integer NOT NULL,
  hgvs_name varchar(255) NOT NULL,
  PRIMARY KEY (variation_id, hgvs_name)
);

--
-- Table: variation_set
--
CREATE TABLE variation_set (
  variation_set_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  name varchar(255),
  description text,
  short_name_attrib_id integer
);

--
-- Table: variation_set_structural_variation
--
CREATE TABLE variation_set_structural_variation (
  structural_variation_id integer NOT NULL,
  variation_set_id integer NOT NULL,
  PRIMARY KEY (structural_variation_id, variation_set_id)
);

--
-- Table: variation_set_structure
--
CREATE TABLE variation_set_structure (
  variation_set_super integer NOT NULL,
  variation_set_sub integer NOT NULL,
  PRIMARY KEY (variation_set_super, variation_set_sub)
);

--
-- Table: variation_set_variation
--
CREATE TABLE variation_set_variation (
  variation_id integer NOT NULL,
  variation_set_id integer NOT NULL,
  PRIMARY KEY (variation_id, variation_set_id)
);

--
-- Table: variation_synonym
--
CREATE TABLE variation_synonym (
  variation_synonym_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  variation_id integer NOT NULL,
  subsnp_id integer,
  source_id integer NOT NULL,
  name varchar(255),
  moltype varchar(50)
);

CREATE UNIQUE INDEX name02 ON variation_synonym (name, source_id);

COMMIT;
