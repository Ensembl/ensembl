
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
);

# Add schema type and schema version to the meta table
INSERT INTO meta (meta_key, meta_value) VALUES 
  ('schema_type', 'ontology'),
  ('schema_version', '73');

# Patches included in this schema file
INSERT INTO meta (meta_key, meta_value)
  VALUES ('patch', 'patch_72_73_a.sql|schema_version');
INSERT INTO meta (meta_key, meta_value)
  VALUES ('patch', 'patch_72_73_b.sql|meta_species');

CREATE TABLE ontology (
  ontology_id   INT UNSIGNED NOT NULL AUTO_INCREMENT,
  name          VARCHAR(64) NOT NULL,
  namespace     VARCHAR(64) NOT NULL,

  PRIMARY KEY (ontology_id),
  UNIQUE INDEX name_namespace_idx (name, namespace)
);

CREATE TABLE subset (
  subset_id     INT UNSIGNED NOT NULL AUTO_INCREMENT,
  name          VARCHAR(64) NOT NULL,
  definition    VARCHAR(128) NOT NULL,

  PRIMARY KEY (subset_id),
  UNIQUE INDEX name_idx (name)
);

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
);

CREATE TABLE synonym (
  synonym_id    INT UNSIGNED NOT NULL AUTO_INCREMENT,
  term_id       INT UNSIGNED NOT NULL,
  name          TEXT NOT NULL,

  PRIMARY KEY (synonym_id),
  UNIQUE INDEX term_synonym_idx (term_id, synonym_id),
  INDEX name_idx (name(50))
);

CREATE TABLE alt_id (
  alt_id        INT UNSIGNED NOT NULL AUTO_INCREMENT,
  term_id       INT UNSIGNED NOT NULL,
  accession     VARCHAR(64) NOT NULL,

  PRIMARY KEY (alt_id),
  UNIQUE INDEX term_alt_idx (term_id, alt_id),
  INDEX accession_idx (accession(50))
);

CREATE TABLE relation_type (
  relation_type_id  INT UNSIGNED NOT NULL AUTO_INCREMENT,
  name              VARCHAR(64) NOT NULL,

  PRIMARY KEY (relation_type_id),
  UNIQUE INDEX name_idx (name)
);

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
);

CREATE TABLE closure (
  closure_id        INT UNSIGNED NOT NULL AUTO_INCREMENT,
  child_term_id     INT UNSIGNED NOT NULL,
  parent_term_id    INT UNSIGNED NOT NULL,
  subparent_term_id INT UNSIGNED,
  distance          TINYINT UNSIGNED NOT NULL,
  ontology_id       INT UNSIGNED NOT NULL,

  PRIMARY KEY (closure_id),
  UNIQUE INDEX child_parent_idx
    (child_term_id, parent_term_id, subparent_term_id, ontology_id),
  INDEX parent_subparent_idx
    (parent_term_id, subparent_term_id)
);

-- There are additional tables in the released databases called
-- "aux_XX_YY_map".  These are created by the "add_subset_maps.pl"
-- scripts.  Please see the README document for further information.

-- $Id$
