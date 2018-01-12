-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2018] EMBL-European Bioinformatics Institute
-- 
-- Licensed under the Apache License, Version 2.0 (the "License");
-- you may not use this file except in compliance with the License.
-- You may obtain a copy of the License at
-- 
--      http://www.apache.org/licenses/LICENSE-2.0
-- 
-- Unless required by applicable law or agreed to in writing, software
-- distributed under the License is distributed on an "AS IS" BASIS,
-- WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
-- See the License for the specific language governing permissions and
-- limitations under the License.


-- ---------------------------------------------------------------------
-- The schema for the ensembl_ontology_NN database.
-- ---------------------------------------------------------------------

CREATE TABLE meta (
  meta_id       INT UNSIGNED NOT NULL AUTO_INCREMENT,
  meta_key      VARCHAR(64) NOT NULL,
  meta_value    VARCHAR(128),
  species_id    INT UNSIGNED DEFAULT NULL,

  PRIMARY KEY (meta_id),
  UNIQUE INDEX key_value_idx (meta_key, meta_value)
) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

# Add schema type and schema version to the meta table
INSERT INTO meta (meta_key, meta_value) VALUES 
  ('schema_type', 'ontology'),
  ('schema_version', '92');

# Patches included in this schema file
INSERT INTO meta (meta_key, meta_value)
  VALUES ('patch', 'patch_91_92_a.sql|schema_version');


CREATE TABLE ontology (
  ontology_id   INT UNSIGNED NOT NULL AUTO_INCREMENT,
  name          VARCHAR(64) NOT NULL,
  namespace     VARCHAR(64) NOT NULL,
  data_version  VARCHAR(64) DEFAULT NULL,

  PRIMARY KEY (ontology_id),
  UNIQUE INDEX name_namespace_idx (name, namespace)
) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

CREATE TABLE subset (
  subset_id     INT UNSIGNED NOT NULL AUTO_INCREMENT,
  name          VARCHAR(64) NOT NULL,
  definition    VARCHAR(128) NOT NULL,

  PRIMARY KEY (subset_id),
  UNIQUE INDEX name_idx (name)
) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

CREATE TABLE term (
  term_id       INT UNSIGNED NOT NULL AUTO_INCREMENT,
  ontology_id   INT UNSIGNED NOT NULL,
  subsets       TEXT,
  accession     VARCHAR(64) NOT NULL,
  name          VARCHAR(255) NOT NULL,
  definition    TEXT,
  is_root       INT NOT NULL DEFAULT 0,
  is_obsolete   INT NOT NULL DEFAULT 0,

  PRIMARY KEY (term_id),
  UNIQUE INDEX accession_idx (accession),
  UNIQUE INDEX ontology_acc_idx (ontology_id, accession),
  INDEX name_idx (name)
) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

CREATE TABLE synonym (
  synonym_id    INT UNSIGNED NOT NULL AUTO_INCREMENT,
  term_id       INT UNSIGNED NOT NULL,
  name          TEXT NOT NULL,
  type		ENUM('EXACT', 'BROAD', 'NARROW', 'RELATED'),
  dbxref        VARCHAR(258)  NULL,

  PRIMARY KEY (synonym_id),
  UNIQUE INDEX term_synonym_idx (term_id, synonym_id),
  INDEX name_idx (name(50))
) COLLATE=utf8_swedish_ci ENGINE=MyISAM;

CREATE TABLE alt_id (
  alt_id        INT UNSIGNED NOT NULL AUTO_INCREMENT,
  term_id       INT UNSIGNED NOT NULL,
  accession     VARCHAR(64) NOT NULL,

  PRIMARY KEY (alt_id),
  UNIQUE INDEX term_alt_idx (term_id, alt_id),
  INDEX accession_idx (accession(50))
) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

CREATE TABLE relation_type (
  relation_type_id  INT UNSIGNED NOT NULL AUTO_INCREMENT,
  name              VARCHAR(64) NOT NULL,

  PRIMARY KEY (relation_type_id),
  UNIQUE INDEX name_idx (name)
) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

CREATE TABLE relation (
  relation_id       INT UNSIGNED NOT NULL AUTO_INCREMENT,
  child_term_id     INT UNSIGNED NOT NULL,
  parent_term_id    INT UNSIGNED NOT NULL,
  relation_type_id  INT UNSIGNED NOT NULL,
  intersection_of   TINYINT UNSIGNED NOT NULL DEFAULT 0,
  ontology_id       INT UNSIGNED NOT NULL,

  PRIMARY KEY (relation_id),
  UNIQUE INDEX child_parent_idx
    (child_term_id, parent_term_id, relation_type_id, intersection_of, ontology_id),
  INDEX parent_idx (parent_term_id)
) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

CREATE TABLE closure (
  closure_id        INT UNSIGNED NOT NULL AUTO_INCREMENT,
  child_term_id     INT UNSIGNED NOT NULL,
  parent_term_id    INT UNSIGNED NOT NULL,
  subparent_term_id INT UNSIGNED,
  distance          TINYINT UNSIGNED NOT NULL,
  ontology_id       INT UNSIGNED NOT NULL,
  confident_relationship BOOL NOT NULL DEFAULT 0,

  PRIMARY KEY (closure_id),
  UNIQUE INDEX child_parent_idx
    (child_term_id, parent_term_id, subparent_term_id, ontology_id),
  INDEX parent_subparent_idx
    (parent_term_id, subparent_term_id)
) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

-- There are additional tables in the released databases called
-- "aux_XX_YY_map".  These are created by the "add_subset_maps.pl"
-- scripts.  Please see the README document for further information.

-- $Id$
