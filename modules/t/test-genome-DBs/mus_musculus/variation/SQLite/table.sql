-- 
-- Created by SQL::Translator::Producer::SQLite
-- Created on Thu Oct 29 15:37:32 2015
-- 

BEGIN TRANSACTION;

--
-- Table: allele
--
CREATE TABLE allele (
  allele_id INTEGER PRIMARY KEY NOT NULL,
  variation_id int(11) NOT NULL,
  subsnp_id int(11) DEFAULT NULL,
  allele_code_id int(11) NOT NULL,
  population_id int(11) DEFAULT NULL,
  frequency float(8,2) DEFAULT NULL,
  count int(11) DEFAULT NULL,
  frequency_submitter_handle int(10) DEFAULT NULL
);

CREATE INDEX variation_idx ON allele (variation_id);

CREATE INDEX subsnp_idx ON allele (subsnp_id);

CREATE INDEX population_idx ON allele (population_id);

--
-- Table: allele_code
--
CREATE TABLE allele_code (
  allele_code_id INTEGER PRIMARY KEY NOT NULL,
  allele varchar(100) DEFAULT NULL
);

CREATE UNIQUE INDEX allele_idx ON allele_code (allele);

--
-- Table: associate_study
--
CREATE TABLE associate_study (
  study1_id int(10) NOT NULL,
  study2_id int(10) NOT NULL,
  PRIMARY KEY (study1_id, study2_id)
);

--
-- Table: attrib
--
CREATE TABLE attrib (
  attrib_id INTEGER PRIMARY KEY NOT NULL DEFAULT 0,
  attrib_type_id smallint(5) NOT NULL DEFAULT 0,
  value text NOT NULL
);

--
-- Table: attrib_set
--
CREATE TABLE attrib_set (
  attrib_set_id int(11) NOT NULL DEFAULT 0,
  attrib_id int(11) NOT NULL DEFAULT 0
);

CREATE INDEX attrib_idx ON attrib_set (attrib_id);

CREATE UNIQUE INDEX set_idx ON attrib_set (attrib_set_id, attrib_id);

--
-- Table: attrib_type
--
CREATE TABLE attrib_type (
  attrib_type_id INTEGER PRIMARY KEY NOT NULL DEFAULT 0,
  code varchar(20) NOT NULL DEFAULT '',
  name varchar(255) NOT NULL DEFAULT '',
  description text
);

CREATE UNIQUE INDEX code_idx ON attrib_type (code);

--
-- Table: compressed_genotype_region
--
CREATE TABLE compressed_genotype_region (
  sample_id int(10) NOT NULL,
  seq_region_id int(10) NOT NULL,
  seq_region_start int(11) NOT NULL,
  seq_region_end int(11) NOT NULL,
  seq_region_strand tinyint(4) NOT NULL,
  genotypes blob
);

CREATE INDEX pos_idx ON compressed_genotype_region (seq_region_id, seq_region_start);

CREATE INDEX sample_idx ON compressed_genotype_region (sample_id);

--
-- Table: compressed_genotype_var
--
CREATE TABLE compressed_genotype_var (
  variation_id int(11) NOT NULL,
  subsnp_id int(11) DEFAULT NULL,
  genotypes blob
);

CREATE INDEX variation_idx02 ON compressed_genotype_var (variation_id);

CREATE INDEX subsnp_idx02 ON compressed_genotype_var (subsnp_id);

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
-- Table: display_group
--
CREATE TABLE display_group (
  display_group_id INTEGER PRIMARY KEY NOT NULL,
  display_priority int(10) NOT NULL,
  display_name varchar(255) NOT NULL
);

CREATE UNIQUE INDEX display_name ON display_group (display_name);

CREATE UNIQUE INDEX display_priority ON display_group (display_priority);

--
-- Table: failed_allele
--
CREATE TABLE failed_allele (
  failed_allele_id INTEGER PRIMARY KEY NOT NULL,
  allele_id int(10) NOT NULL,
  failed_description_id int(10) NOT NULL
);

CREATE UNIQUE INDEX allele_idx02 ON failed_allele (allele_id, failed_description_id);

--
-- Table: failed_description
--
CREATE TABLE failed_description (
  failed_description_id INTEGER PRIMARY KEY NOT NULL,
  description text NOT NULL
);

--
-- Table: failed_structural_variation
--
CREATE TABLE failed_structural_variation (
  failed_structural_variation_id INTEGER PRIMARY KEY NOT NULL,
  structural_variation_id int(10) NOT NULL,
  failed_description_id int(10) NOT NULL
);

CREATE UNIQUE INDEX structural_variation_idx ON failed_structural_variation (structural_variation_id, failed_description_id);

--
-- Table: failed_variation
--
CREATE TABLE failed_variation (
  failed_variation_id INTEGER PRIMARY KEY NOT NULL,
  variation_id int(10) NOT NULL,
  failed_description_id int(10) NOT NULL
);

CREATE UNIQUE INDEX variation_idx03 ON failed_variation (variation_id, failed_description_id);

--
-- Table: genotype_code
--
CREATE TABLE genotype_code (
  genotype_code_id int(11) NOT NULL,
  allele_code_id int(11) NOT NULL,
  haplotype_id tinyint(2) NOT NULL,
  phased tinyint(2) DEFAULT NULL
);

CREATE INDEX genotype_code_id ON genotype_code (genotype_code_id);

CREATE INDEX allele_code_id ON genotype_code (allele_code_id);

--
-- Table: individual
--
CREATE TABLE individual (
  individual_id INTEGER PRIMARY KEY NOT NULL,
  name varchar(255) DEFAULT NULL,
  description text,
  gender enum(7) NOT NULL DEFAULT 'Unknown',
  father_individual_id int(10) DEFAULT NULL,
  mother_individual_id int(10) DEFAULT NULL,
  individual_type_id int(10) NOT NULL DEFAULT 0
);

CREATE INDEX father_individual_idx ON individual (father_individual_id);

CREATE INDEX mother_individual_idx ON individual (mother_individual_id);

--
-- Table: individual_synonym
--
CREATE TABLE individual_synonym (
  synonym_id INTEGER PRIMARY KEY NOT NULL,
  individual_id int(10) NOT NULL,
  source_id int(10) NOT NULL,
  name varchar(255) DEFAULT NULL
);

CREATE INDEX individual_idx ON individual_synonym (individual_id);

CREATE INDEX name ON individual_synonym (name, source_id);

--
-- Table: individual_type
--
CREATE TABLE individual_type (
  individual_type_id INTEGER PRIMARY KEY NOT NULL,
  name varchar(255) NOT NULL,
  description text
);

--
-- Table: meta
--
CREATE TABLE meta (
  meta_id INTEGER PRIMARY KEY NOT NULL,
  species_id int(10) DEFAULT 1,
  meta_key varchar(40) NOT NULL,
  meta_value varchar(255) NOT NULL
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

CREATE UNIQUE INDEX table_name ON meta_coord (table_name, coord_system_id);

--
-- Table: motif_feature_variation
--
CREATE TABLE motif_feature_variation (
  motif_feature_variation_id INTEGER PRIMARY KEY NOT NULL,
  variation_feature_id int(11) NOT NULL,
  feature_stable_id varchar(128) DEFAULT NULL,
  motif_feature_id int(11) NOT NULL,
  allele_string text,
  somatic tinyint(1) NOT NULL DEFAULT 0,
  consequence_types varchar(23) DEFAULT NULL,
  motif_name varchar(60) DEFAULT NULL,
  motif_start int(11) DEFAULT NULL,
  motif_end int(11) DEFAULT NULL,
  motif_score_delta float(8,2) DEFAULT NULL,
  in_informative_position tinyint(1) NOT NULL DEFAULT 0
);

CREATE INDEX variation_feature_idx ON motif_feature_variation (variation_feature_id);

CREATE INDEX consequence_type_idx ON motif_feature_variation (consequence_types);

CREATE INDEX somatic_feature_idx ON motif_feature_variation (feature_stable_id, somatic);

--
-- Table: phenotype
--
CREATE TABLE phenotype (
  phenotype_id INTEGER PRIMARY KEY NOT NULL,
  stable_id varchar(255) DEFAULT NULL,
  name varchar(50) DEFAULT NULL,
  description varchar(255) DEFAULT NULL
);

CREATE INDEX name_idx02 ON phenotype (name);

CREATE INDEX stable_idx ON phenotype (stable_id);

CREATE UNIQUE INDEX desc_idx ON phenotype (description);

--
-- Table: phenotype_feature
--
CREATE TABLE phenotype_feature (
  phenotype_feature_id INTEGER PRIMARY KEY NOT NULL,
  phenotype_id int(11) DEFAULT NULL,
  source_id int(11) DEFAULT NULL,
  study_id int(11) DEFAULT NULL,
  type enum(29) DEFAULT NULL,
  object_id varchar(255) DEFAULT NULL,
  is_significant tinyint(1) DEFAULT 1,
  seq_region_id int(11) DEFAULT NULL,
  seq_region_start int(11) DEFAULT NULL,
  seq_region_end int(11) DEFAULT NULL,
  seq_region_strand tinyint(4) DEFAULT NULL
);

CREATE INDEX phenotype_idx ON phenotype_feature (phenotype_id);

CREATE INDEX object_idx ON phenotype_feature (object_id, type);

CREATE INDEX type_idx ON phenotype_feature (type);

CREATE INDEX pos_idx02 ON phenotype_feature (seq_region_id, seq_region_start, seq_region_end);

CREATE INDEX source_idx ON phenotype_feature (source_id);

--
-- Table: phenotype_feature_attrib
--
CREATE TABLE phenotype_feature_attrib (
  phenotype_feature_id int(11) NOT NULL,
  attrib_type_id int(11) DEFAULT NULL,
  value varchar(255) DEFAULT NULL
);

CREATE INDEX phenotype_feature_idx ON phenotype_feature_attrib (phenotype_feature_id);

CREATE INDEX type_value_idx ON phenotype_feature_attrib (attrib_type_id, value);

--
-- Table: population
--
CREATE TABLE population (
  population_id INTEGER PRIMARY KEY NOT NULL,
  name varchar(255) DEFAULT NULL,
  size int(10) DEFAULT NULL,
  description text,
  collection tinyint(1) DEFAULT 0,
  freqs_from_gts tinyint(1) DEFAULT NULL,
  display enum(15) DEFAULT 'UNDISPLAYABLE',
  display_group_id tinyint(1) DEFAULT NULL
);

CREATE INDEX name_idx03 ON population (name);

--
-- Table: population_genotype
--
CREATE TABLE population_genotype (
  population_genotype_id INTEGER PRIMARY KEY NOT NULL,
  variation_id int(11) NOT NULL,
  subsnp_id int(11) DEFAULT NULL,
  genotype_code_id int(11) DEFAULT NULL,
  frequency float(8,2) DEFAULT NULL,
  population_id int(10) DEFAULT NULL,
  count int(10) DEFAULT NULL
);

CREATE INDEX population_idx02 ON population_genotype (population_id);

CREATE INDEX variation_idx04 ON population_genotype (variation_id);

CREATE INDEX subsnp_idx03 ON population_genotype (subsnp_id);

--
-- Table: population_structure
--
CREATE TABLE population_structure (
  super_population_id int(10) NOT NULL,
  sub_population_id int(10) NOT NULL
);

CREATE INDEX sub_population_idx ON population_structure (sub_population_id);

CREATE UNIQUE INDEX super_population_idx ON population_structure (super_population_id, sub_population_id);

--
-- Table: population_synonym
--
CREATE TABLE population_synonym (
  synonym_id INTEGER PRIMARY KEY NOT NULL,
  population_id int(10) NOT NULL,
  source_id int(10) NOT NULL,
  name varchar(255) DEFAULT NULL
);

CREATE INDEX population_idx03 ON population_synonym (population_id);

CREATE INDEX name02 ON population_synonym (name, source_id);

--
-- Table: protein_function_predictions
--
CREATE TABLE protein_function_predictions (
  translation_md5_id int(11) NOT NULL,
  analysis_attrib_id int(11) NOT NULL,
  prediction_matrix mediumblob(16777215),
  PRIMARY KEY (translation_md5_id, analysis_attrib_id)
);

--
-- Table: protein_function_predictions_attrib
--
CREATE TABLE protein_function_predictions_attrib (
  translation_md5_id int(11) NOT NULL,
  analysis_attrib_id int(11) NOT NULL,
  attrib_type_id int(11) NOT NULL,
  position_values blob,
  PRIMARY KEY (translation_md5_id, analysis_attrib_id, attrib_type_id)
);

--
-- Table: publication
--
CREATE TABLE publication (
  publication_id INTEGER PRIMARY KEY NOT NULL,
  title varchar(255) DEFAULT NULL,
  authors varchar(255) DEFAULT NULL,
  pmid int(10) DEFAULT NULL,
  pmcid varchar(255) DEFAULT NULL,
  year int(10) DEFAULT NULL,
  doi varchar(50) DEFAULT NULL,
  ucsc_id varchar(50) DEFAULT NULL
);

CREATE INDEX pmid_idx ON publication (pmid);

CREATE INDEX doi_idx ON publication (doi);

--
-- Table: read_coverage
--
CREATE TABLE read_coverage (
  seq_region_id int(10) NOT NULL,
  seq_region_start int(11) NOT NULL,
  seq_region_end int(11) NOT NULL,
  level tinyint(4) NOT NULL,
  sample_id int(10) NOT NULL
);

CREATE INDEX seq_region_idx ON read_coverage (seq_region_id, seq_region_start);

CREATE INDEX sample_idx02 ON read_coverage (sample_id);

--
-- Table: regulatory_feature_variation
--
CREATE TABLE regulatory_feature_variation (
  regulatory_feature_variation_id INTEGER PRIMARY KEY NOT NULL,
  variation_feature_id int(11) NOT NULL,
  feature_stable_id varchar(128) DEFAULT NULL,
  feature_type text,
  allele_string text,
  somatic tinyint(1) NOT NULL DEFAULT 0,
  consequence_types varchar(31) DEFAULT NULL
);

CREATE INDEX variation_feature_idx02 ON regulatory_feature_variation (variation_feature_id);

CREATE INDEX consequence_type_idx02 ON regulatory_feature_variation (consequence_types);

CREATE INDEX somatic_feature_idx02 ON regulatory_feature_variation (feature_stable_id, somatic);

--
-- Table: sample
--
CREATE TABLE sample (
  sample_id INTEGER PRIMARY KEY NOT NULL,
  individual_id int(10) NOT NULL,
  name varchar(255) DEFAULT NULL,
  description text,
  study_id int(10) DEFAULT NULL,
  display enum(15) DEFAULT 'UNDISPLAYABLE',
  has_coverage tinyint(1) NOT NULL DEFAULT 0,
  variation_set_id varchar(2) DEFAULT NULL
);

CREATE INDEX individual_idx02 ON sample (individual_id);

CREATE INDEX study_idx ON sample (study_id);

--
-- Table: sample_genotype_multiple_bp
--
CREATE TABLE sample_genotype_multiple_bp (
  variation_id int(10) NOT NULL,
  subsnp_id int(15) DEFAULT NULL,
  allele_1 varchar(100) DEFAULT NULL,
  allele_2 varchar(100) DEFAULT NULL,
  sample_id int(10) DEFAULT NULL
);

CREATE INDEX variation_idx05 ON sample_genotype_multiple_bp (variation_id);

CREATE INDEX subsnp_idx04 ON sample_genotype_multiple_bp (subsnp_id);

CREATE INDEX sample_idx03 ON sample_genotype_multiple_bp (sample_id);

--
-- Table: sample_population
--
CREATE TABLE sample_population (
  sample_id int(10) NOT NULL,
  population_id int(10) NOT NULL
);

CREATE INDEX population_idx04 ON sample_population (population_id);

CREATE INDEX sample_idx04 ON sample_population (sample_id);

--
-- Table: seq_region
--
CREATE TABLE seq_region (
  seq_region_id INTEGER PRIMARY KEY NOT NULL,
  name varchar(40) NOT NULL,
  coord_system_id int(10) NOT NULL
);

CREATE INDEX cs_idx ON seq_region (coord_system_id);

CREATE UNIQUE INDEX name_cs_idx ON seq_region (name, coord_system_id);

--
-- Table: source
--
CREATE TABLE source (
  source_id INTEGER PRIMARY KEY NOT NULL,
  name varchar(24) NOT NULL,
  version int(11) DEFAULT NULL,
  description varchar(255) DEFAULT NULL,
  url varchar(255) DEFAULT NULL,
  type enum(4) DEFAULT NULL,
  somatic_status enum(8) DEFAULT 'germline',
  data_types varchar(20) DEFAULT NULL
);

--
-- Table: structural_variation
--
CREATE TABLE structural_variation (
  structural_variation_id INTEGER PRIMARY KEY NOT NULL,
  variation_name varchar(255) DEFAULT NULL,
  alias varchar(255) DEFAULT NULL,
  source_id int(10) NOT NULL,
  study_id int(10) DEFAULT NULL,
  class_attrib_id int(10) NOT NULL DEFAULT 0,
  clinical_significance varchar(22) DEFAULT NULL,
  validation_status enum(13) DEFAULT NULL,
  is_evidence tinyint(4) DEFAULT 0,
  somatic tinyint(1) NOT NULL DEFAULT 0,
  copy_number tinyint(2) DEFAULT NULL
);

CREATE INDEX study_idx02 ON structural_variation (study_id);

CREATE INDEX attrib_idx02 ON structural_variation (class_attrib_id);

CREATE INDEX source_idx02 ON structural_variation (source_id);

CREATE UNIQUE INDEX variation_name ON structural_variation (variation_name);

--
-- Table: structural_variation_association
--
CREATE TABLE structural_variation_association (
  structural_variation_id int(10) NOT NULL,
  supporting_structural_variation_id int(10) NOT NULL,
  PRIMARY KEY (structural_variation_id, supporting_structural_variation_id)
);

CREATE INDEX structural_variation_idx02 ON structural_variation_association (structural_variation_id);

CREATE INDEX supporting_structural_variation_idx ON structural_variation_association (supporting_structural_variation_id);

--
-- Table: structural_variation_feature
--
CREATE TABLE structural_variation_feature (
  structural_variation_feature_id INTEGER PRIMARY KEY NOT NULL,
  seq_region_id int(10) NOT NULL,
  outer_start int(11) DEFAULT NULL,
  seq_region_start int(11) NOT NULL,
  inner_start int(11) DEFAULT NULL,
  inner_end int(11) DEFAULT NULL,
  seq_region_end int(11) NOT NULL,
  outer_end int(11) DEFAULT NULL,
  seq_region_strand tinyint(4) NOT NULL,
  structural_variation_id int(10) NOT NULL,
  variation_name varchar(255) DEFAULT NULL,
  source_id int(10) NOT NULL,
  study_id int(10) DEFAULT NULL,
  class_attrib_id int(10) NOT NULL DEFAULT 0,
  allele_string longtext(4294967295),
  is_evidence tinyint(1) NOT NULL DEFAULT 0,
  variation_set_id varchar(2) NOT NULL DEFAULT '',
  somatic tinyint(1) NOT NULL DEFAULT 0,
  breakpoint_order tinyint(4) DEFAULT NULL,
  length int(10) DEFAULT NULL
);

CREATE INDEX pos_idx03 ON structural_variation_feature (seq_region_id, seq_region_start, seq_region_end);

CREATE INDEX structural_variation_idx03 ON structural_variation_feature (structural_variation_id);

CREATE INDEX attrib_idx03 ON structural_variation_feature (class_attrib_id);

CREATE INDEX source_idx03 ON structural_variation_feature (source_id);

CREATE INDEX variation_set_idx ON structural_variation_feature (variation_set_id);

CREATE INDEX study_idx03 ON structural_variation_feature (study_id);

--
-- Table: structural_variation_sample
--
CREATE TABLE structural_variation_sample (
  structural_variation_sample_id INTEGER PRIMARY KEY NOT NULL,
  structural_variation_id int(10) NOT NULL,
  sample_id int(10) DEFAULT NULL
);

CREATE INDEX structural_variation_idx04 ON structural_variation_sample (structural_variation_id);

CREATE INDEX sample_idx05 ON structural_variation_sample (sample_id);

--
-- Table: study
--
CREATE TABLE study (
  study_id INTEGER PRIMARY KEY NOT NULL,
  source_id int(10) NOT NULL,
  name varchar(255) DEFAULT NULL,
  description text,
  url varchar(255) DEFAULT NULL,
  external_reference varchar(255) DEFAULT NULL,
  study_type varchar(255) DEFAULT NULL
);

CREATE INDEX source_idx04 ON study (source_id);

--
-- Table: submitter_handle
--
CREATE TABLE submitter_handle (
  handle_id INTEGER PRIMARY KEY NOT NULL,
  handle varchar(25) DEFAULT NULL
);

CREATE UNIQUE INDEX handle ON submitter_handle (handle);

--
-- Table: subsnp_handle
--
CREATE TABLE subsnp_handle (
  subsnp_id INTEGER PRIMARY KEY NOT NULL,
  handle varchar(20) DEFAULT NULL
);

--
-- Table: subsnp_map
--
CREATE TABLE subsnp_map (
  variation_id int(11) NOT NULL,
  subsnp_id int(11) DEFAULT NULL
);

CREATE INDEX variation_idx06 ON subsnp_map (variation_id);

--
-- Table: tagged_variation_feature
--
CREATE TABLE tagged_variation_feature (
  variation_feature_id int(10) NOT NULL,
  tagged_variation_feature_id int(10) DEFAULT NULL,
  population_id int(10) NOT NULL
);

CREATE INDEX tag_idx ON tagged_variation_feature (variation_feature_id);

CREATE INDEX tagged_idx ON tagged_variation_feature (tagged_variation_feature_id);

CREATE INDEX population_idx05 ON tagged_variation_feature (population_id);

--
-- Table: tmp_sample_genotype_single_bp
--
CREATE TABLE tmp_sample_genotype_single_bp (
  variation_id int(10) NOT NULL,
  subsnp_id int(15) DEFAULT NULL,
  allele_1 char(1) DEFAULT NULL,
  allele_2 char(1) DEFAULT NULL,
  sample_id int(10) NOT NULL
);

CREATE INDEX variation_idx07 ON tmp_sample_genotype_single_bp (variation_id);

CREATE INDEX subsnp_idx05 ON tmp_sample_genotype_single_bp (subsnp_id);

CREATE INDEX sample_idx06 ON tmp_sample_genotype_single_bp (sample_id);

--
-- Table: transcript_variation
--
CREATE TABLE transcript_variation (
  transcript_variation_id INTEGER PRIMARY KEY NOT NULL,
  variation_feature_id int(11) NOT NULL,
  feature_stable_id varchar(128) DEFAULT NULL,
  allele_string text,
  somatic tinyint(1) NOT NULL DEFAULT 0,
  consequence_types varchar(34) DEFAULT NULL,
  cds_start int(11) DEFAULT NULL,
  cds_end int(11) DEFAULT NULL,
  cdna_start int(11) DEFAULT NULL,
  cdna_end int(11) DEFAULT NULL,
  translation_start int(11) DEFAULT NULL,
  translation_end int(11) DEFAULT NULL,
  distance_to_transcript int(11) DEFAULT NULL,
  codon_allele_string text,
  pep_allele_string text,
  hgvs_genomic text,
  hgvs_transcript text,
  hgvs_protein text,
  polyphen_prediction enum(17) DEFAULT NULL,
  sift_prediction enum(28) DEFAULT NULL,
  polyphen_score float(8,2) DEFAULT NULL,
  sift_score float(8,2) DEFAULT NULL,
  display int(1) DEFAULT 1
);

CREATE INDEX variation_feature_idx03 ON transcript_variation (variation_feature_id);

CREATE INDEX consequence_type_idx03 ON transcript_variation (consequence_types);

CREATE INDEX somatic_feature_idx03 ON transcript_variation (feature_stable_id, somatic);

--
-- Table: translation_md5
--
CREATE TABLE translation_md5 (
  translation_md5_id INTEGER PRIMARY KEY NOT NULL,
  translation_md5 char(32) NOT NULL
);

CREATE UNIQUE INDEX md5_idx ON translation_md5 (translation_md5);

--
-- Table: variation
--
CREATE TABLE variation (
  variation_id INTEGER PRIMARY KEY NOT NULL,
  source_id int(10) NOT NULL,
  name varchar(255) DEFAULT NULL,
  validation_status varchar(10) DEFAULT NULL,
  ancestral_allele varchar(255) DEFAULT NULL,
  flipped tinyint(1) DEFAULT NULL,
  class_attrib_id int(10) DEFAULT 0,
  somatic tinyint(1) NOT NULL DEFAULT 0,
  minor_allele varchar(50) DEFAULT NULL,
  minor_allele_freq float(8,2) DEFAULT NULL,
  minor_allele_count int(10) DEFAULT NULL,
  clinical_significance varchar(22) DEFAULT NULL,
  evidence_attribs varchar(3) DEFAULT NULL,
  display int(1) DEFAULT 1
);

CREATE INDEX source_idx05 ON variation (source_id);

CREATE UNIQUE INDEX name03 ON variation (name);

--
-- Table: variation_attrib
--
CREATE TABLE variation_attrib (
  variation_id int(11) NOT NULL,
  attrib_id int(11) DEFAULT NULL,
  value varchar(255) DEFAULT NULL
);

CREATE INDEX variation_idx08 ON variation_attrib (variation_id);

CREATE INDEX attrib_value_idx ON variation_attrib (attrib_id, value);

--
-- Table: variation_citation
--
CREATE TABLE variation_citation (
  variation_id int(10) NOT NULL,
  publication_id int(10) NOT NULL,
  PRIMARY KEY (variation_id, publication_id)
);

--
-- Table: variation_feature
--
CREATE TABLE variation_feature (
  variation_feature_id INTEGER PRIMARY KEY NOT NULL,
  seq_region_id int(10) NOT NULL,
  seq_region_start int(11) NOT NULL,
  seq_region_end int(11) NOT NULL,
  seq_region_strand tinyint(4) NOT NULL,
  variation_id int(10) NOT NULL,
  allele_string varchar(100) DEFAULT NULL,
  variation_name varchar(255) DEFAULT NULL,
  map_weight int(11) NOT NULL,
  flags varchar(9) DEFAULT NULL,
  source_id int(10) NOT NULL,
  validation_status varchar(10) DEFAULT NULL,
  consequence_types varchar(34) NOT NULL DEFAULT 'intergenic_variant',
  variation_set_id varchar(2) NOT NULL DEFAULT '',
  class_attrib_id int(10) DEFAULT 0,
  somatic tinyint(1) NOT NULL DEFAULT 0,
  minor_allele varchar(50) DEFAULT NULL,
  minor_allele_freq float(8,2) DEFAULT NULL,
  minor_allele_count int(10) DEFAULT NULL,
  alignment_quality double(8,2) DEFAULT NULL,
  evidence_attribs varchar(3) DEFAULT NULL,
  clinical_significance varchar(22) DEFAULT NULL,
  display int(1) DEFAULT 1
);

CREATE INDEX pos_idx04 ON variation_feature (seq_region_id, seq_region_start, seq_region_end);

CREATE INDEX variation_idx09 ON variation_feature (variation_id);

CREATE INDEX variation_set_idx02 ON variation_feature (variation_set_id);

CREATE INDEX consequence_type_idx04 ON variation_feature (consequence_types);

CREATE INDEX source_idx06 ON variation_feature (source_id);

--
-- Table: variation_genename
--
CREATE TABLE variation_genename (
  variation_id int(10) NOT NULL,
  gene_name varchar(255) NOT NULL,
  PRIMARY KEY (variation_id, gene_name)
);

--
-- Table: variation_hgvs
--
CREATE TABLE variation_hgvs (
  variation_id int(10) NOT NULL,
  hgvs_name varchar(255) NOT NULL,
  PRIMARY KEY (variation_id, hgvs_name)
);

--
-- Table: variation_set
--
CREATE TABLE variation_set (
  variation_set_id INTEGER PRIMARY KEY NOT NULL,
  name varchar(255) DEFAULT NULL,
  description text,
  short_name_attrib_id int(10) DEFAULT NULL
);

CREATE INDEX name_idx04 ON variation_set (name);

--
-- Table: variation_set_structural_variation
--
CREATE TABLE variation_set_structural_variation (
  structural_variation_id int(10) NOT NULL,
  variation_set_id int(10) NOT NULL,
  PRIMARY KEY (structural_variation_id, variation_set_id)
);

--
-- Table: variation_set_structure
--
CREATE TABLE variation_set_structure (
  variation_set_super int(10) NOT NULL,
  variation_set_sub int(10) NOT NULL,
  PRIMARY KEY (variation_set_super, variation_set_sub)
);

CREATE INDEX sub_idx ON variation_set_structure (variation_set_sub, variation_set_super);

--
-- Table: variation_set_variation
--
CREATE TABLE variation_set_variation (
  variation_id int(10) NOT NULL,
  variation_set_id int(10) NOT NULL,
  PRIMARY KEY (variation_id, variation_set_id)
);

CREATE INDEX variation_set_idx03 ON variation_set_variation (variation_set_id, variation_id);

--
-- Table: variation_synonym
--
CREATE TABLE variation_synonym (
  variation_synonym_id INTEGER PRIMARY KEY NOT NULL,
  variation_id int(10) NOT NULL,
  subsnp_id int(15) DEFAULT NULL,
  source_id int(10) NOT NULL,
  name varchar(255) DEFAULT NULL,
  moltype varchar(50) DEFAULT NULL
);

CREATE INDEX variation_idx10 ON variation_synonym (variation_id);

CREATE INDEX subsnp_idx06 ON variation_synonym (subsnp_id);

CREATE INDEX source_idx07 ON variation_synonym (source_id);

CREATE UNIQUE INDEX name04 ON variation_synonym (name, source_id);

COMMIT;
