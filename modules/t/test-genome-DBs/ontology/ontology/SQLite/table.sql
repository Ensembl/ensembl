-- 
-- Created by SQL::Translator::Producer::SQLite
-- Created on Fri Jan 12 13:38:37 2018
-- 

BEGIN TRANSACTION;

--
-- Table: alt_id
--
CREATE TABLE alt_id (
  alt_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  term_id integer NOT NULL,
  accession varchar(64) NOT NULL
);

CREATE UNIQUE INDEX term_alt_idx ON alt_id (term_id, alt_id);

--
-- Table: aux_GO_Cross_product_review_map
--
CREATE TABLE aux_GO_Cross_product_review_map (
  term_id integer NOT NULL,
  subset_term_id integer NOT NULL,
  distance tinyint NOT NULL
);

CREATE UNIQUE INDEX map_idx ON aux_GO_Cross_product_review_map (term_id, subset_term_id);

--
-- Table: aux_GO_goslim_aspergillus_map
--
CREATE TABLE aux_GO_goslim_aspergillus_map (
  term_id integer NOT NULL,
  subset_term_id integer NOT NULL,
  distance tinyint NOT NULL
);

CREATE UNIQUE INDEX map_idx02 ON aux_GO_goslim_aspergillus_map (term_id, subset_term_id);

--
-- Table: aux_GO_goslim_candida_map
--
CREATE TABLE aux_GO_goslim_candida_map (
  term_id integer NOT NULL,
  subset_term_id integer NOT NULL,
  distance tinyint NOT NULL
);

CREATE UNIQUE INDEX map_idx03 ON aux_GO_goslim_candida_map (term_id, subset_term_id);

--
-- Table: aux_GO_goslim_generic_map
--
CREATE TABLE aux_GO_goslim_generic_map (
  term_id integer NOT NULL,
  subset_term_id integer NOT NULL,
  distance tinyint NOT NULL
);

CREATE UNIQUE INDEX map_idx04 ON aux_GO_goslim_generic_map (term_id, subset_term_id);

--
-- Table: aux_GO_goslim_metagenomics_map
--
CREATE TABLE aux_GO_goslim_metagenomics_map (
  term_id integer NOT NULL,
  subset_term_id integer NOT NULL,
  distance tinyint NOT NULL
);

CREATE UNIQUE INDEX map_idx05 ON aux_GO_goslim_metagenomics_map (term_id, subset_term_id);

--
-- Table: aux_GO_goslim_pir_map
--
CREATE TABLE aux_GO_goslim_pir_map (
  term_id integer NOT NULL,
  subset_term_id integer NOT NULL,
  distance tinyint NOT NULL
);

CREATE UNIQUE INDEX map_idx06 ON aux_GO_goslim_pir_map (term_id, subset_term_id);

--
-- Table: aux_GO_goslim_plant_map
--
CREATE TABLE aux_GO_goslim_plant_map (
  term_id integer NOT NULL,
  subset_term_id integer NOT NULL,
  distance tinyint NOT NULL
);

CREATE UNIQUE INDEX map_idx07 ON aux_GO_goslim_plant_map (term_id, subset_term_id);

--
-- Table: aux_GO_goslim_pombe_map
--
CREATE TABLE aux_GO_goslim_pombe_map (
  term_id integer NOT NULL,
  subset_term_id integer NOT NULL,
  distance tinyint NOT NULL
);

CREATE UNIQUE INDEX map_idx08 ON aux_GO_goslim_pombe_map (term_id, subset_term_id);

--
-- Table: aux_GO_goslim_yeast_map
--
CREATE TABLE aux_GO_goslim_yeast_map (
  term_id integer NOT NULL,
  subset_term_id integer NOT NULL,
  distance tinyint NOT NULL
);

CREATE UNIQUE INDEX map_idx09 ON aux_GO_goslim_yeast_map (term_id, subset_term_id);

--
-- Table: aux_GO_gosubset_prok_map
--
CREATE TABLE aux_GO_gosubset_prok_map (
  term_id integer NOT NULL,
  subset_term_id integer NOT NULL,
  distance tinyint NOT NULL
);

CREATE UNIQUE INDEX map_idx10 ON aux_GO_gosubset_prok_map (term_id, subset_term_id);

--
-- Table: aux_GO_high_level_annotation_qc_map
--
CREATE TABLE aux_GO_high_level_annotation_qc_map (
  term_id integer NOT NULL,
  subset_term_id integer NOT NULL,
  distance tinyint NOT NULL
);

CREATE UNIQUE INDEX map_idx11 ON aux_GO_high_level_annotation_qc_map (term_id, subset_term_id);

--
-- Table: aux_GO_mf_needs_review_map
--
CREATE TABLE aux_GO_mf_needs_review_map (
  term_id integer NOT NULL,
  subset_term_id integer NOT NULL,
  distance tinyint NOT NULL
);

CREATE UNIQUE INDEX map_idx12 ON aux_GO_mf_needs_review_map (term_id, subset_term_id);

--
-- Table: aux_GO_virus_checked_map
--
CREATE TABLE aux_GO_virus_checked_map (
  term_id integer NOT NULL,
  subset_term_id integer NOT NULL,
  distance tinyint NOT NULL
);

CREATE UNIQUE INDEX map_idx13 ON aux_GO_virus_checked_map (term_id, subset_term_id);

--
-- Table: aux_SO_DBVAR_map
--
CREATE TABLE aux_SO_DBVAR_map (
  term_id integer NOT NULL,
  subset_term_id integer NOT NULL,
  distance tinyint NOT NULL
);

CREATE UNIQUE INDEX map_idx14 ON aux_SO_DBVAR_map (term_id, subset_term_id);

--
-- Table: aux_SO_SOFA_map
--
CREATE TABLE aux_SO_SOFA_map (
  term_id integer NOT NULL,
  subset_term_id integer NOT NULL,
  distance tinyint NOT NULL
);

CREATE UNIQUE INDEX map_idx15 ON aux_SO_SOFA_map (term_id, subset_term_id);

--
-- Table: aux_SO_biosapiens_map
--
CREATE TABLE aux_SO_biosapiens_map (
  term_id integer NOT NULL,
  subset_term_id integer NOT NULL,
  distance tinyint NOT NULL
);

CREATE UNIQUE INDEX map_idx16 ON aux_SO_biosapiens_map (term_id, subset_term_id);

--
-- Table: closure
--
CREATE TABLE closure (
  closure_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  child_term_id integer NOT NULL,
  parent_term_id integer NOT NULL,
  subparent_term_id integer,
  distance tinyint NOT NULL,
  ontology_id integer NOT NULL,
  confident_relationship tinyint NOT NULL DEFAULT 0
);

CREATE UNIQUE INDEX child_parent_idx ON closure (child_term_id, parent_term_id, subparent_term_id, ontology_id);

--
-- Table: meta
--
CREATE TABLE meta (
  meta_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  meta_key varchar(64) NOT NULL,
  meta_value varchar(128),
  species_id integer
);

CREATE UNIQUE INDEX key_value_idx ON meta (meta_key, meta_value);

--
-- Table: ontology
--
CREATE TABLE ontology (
  ontology_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  name varchar(64) NOT NULL,
  namespace varchar(64) NOT NULL,
  data_version varchar(64)
);

CREATE UNIQUE INDEX name_namespace_idx ON ontology (name, namespace);

--
-- Table: relation
--
CREATE TABLE relation (
  relation_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  child_term_id integer NOT NULL,
  parent_term_id integer NOT NULL,
  relation_type_id integer NOT NULL,
  intersection_of tinyint NOT NULL DEFAULT 0,
  ontology_id integer NOT NULL
);

CREATE UNIQUE INDEX child_parent_idx02 ON relation (child_term_id, parent_term_id, relation_type_id, intersection_of, ontology_id);

--
-- Table: relation_type
--
CREATE TABLE relation_type (
  relation_type_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  name varchar(64) NOT NULL
);

CREATE UNIQUE INDEX name_idx ON relation_type (name);

--
-- Table: subset
--
CREATE TABLE subset (
  subset_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  name varchar(64) NOT NULL,
  definition varchar(128) NOT NULL
);

CREATE UNIQUE INDEX name_idx02 ON subset (name);

--
-- Table: synonym
--
CREATE TABLE synonym (
  synonym_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  term_id integer NOT NULL,
  name mediumtext NOT NULL,
  type enum,
  dbxref varchar(256) NOT NULL
);

CREATE UNIQUE INDEX term_synonym_idx ON synonym (term_id, synonym_id);

--
-- Table: term
--
CREATE TABLE term (
  term_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  ontology_id integer NOT NULL,
  subsets text,
  accession varchar(64) NOT NULL,
  name varchar(255) NOT NULL,
  definition text,
  is_root integer,
  is_obsolete integer
);

CREATE UNIQUE INDEX accession_idx ON term (accession);

CREATE UNIQUE INDEX ontology_acc_idx ON term (ontology_id, accession);

COMMIT;
