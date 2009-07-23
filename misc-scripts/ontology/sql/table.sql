-- ---------------------------------------------------------------------
-- The schema for the ensembl_ontology_NN database.
-- ---------------------------------------------------------------------

DROP TABLE IF EXISTS meta;
CREATE TABLE meta (
  meta_id       INT UNSIGNED NOT NULL AUTO_INCREMENT,
  meta_key      VARCHAR(64) NOT NULL,
  meta_value    VARCHAR(128),

  PRIMARY KEY (meta_id),
  UNIQUE INDEX key_value_idx (meta_key, meta_value)
);

DROP TABLE IF EXISTS ontology;
CREATE TABLE ontology (
  ontology_id   INT UNSIGNED NOT NULL AUTO_INCREMENT,
  name          VARCHAR(64) NOT NULL,
  namespace     VARCHAR(64) NOT NULL,

  PRIMARY KEY (ontology_id),
  UNIQUE INDEX name_namespace_idx (name, namespace)
);

DROP TABLE IF EXISTS subset;
CREATE TABLE subset (
  subset_id     INT UNSIGNED NOT NULL AUTO_INCREMENT,
  name          VARCHAR(64) NOT NULL,
  definition    VARCHAR(128) NOT NULL,

  PRIMARY KEY (subset_id),
  UNIQUE INDEX name_idx (name)
);

DROP TABLE IF EXISTS term;
CREATE TABLE term (
  term_id       INT UNSIGNED NOT NULL AUTO_INCREMENT,
  ontology_id   INT UNSIGNED NOT NULL,
  subsets       SET(
                  'SOFA', 'biosapiens', -- SO
                  'gosubset_prok', 'goslim_candida',               -- GO
                  'goslim_generic', 'goslim_goa', 'goslim_pir',    -- GO
                  'goslim_plant', 'goslim_pombe', 'goslim_yeast'   -- GO
                ),
  accession     CHAR(10) NOT NULL,
  name          VARCHAR(255) NOT NULL,
  definition    TEXT,

  PRIMARY KEY (term_id),
  INDEX ontology_idx (ontology_id),
  UNIQUE INDEX acc_ontology_idx (accession, ontology_id)
);

DROP TABLE IF EXISTS relation_type;
CREATE TABLE relation_type (
  relation_type_id  INT UNSIGNED NOT NULL AUTO_INCREMENT,
  name              VARCHAR(64) NOT NULL,

  PRIMARY KEY (relation_type_id),
  UNIQUE INDEX name_idx (name)
);

DROP TABLE IF EXISTS relation;
CREATE TABLE relation (
  relation_id       INT UNSIGNED NOT NULL AUTO_INCREMENT,
  child_term_id     INT UNSIGNED NOT NULL,
  parent_term_id    INT UNSIGNED NOT NULL,
  relation_type_id  INT UNSIGNED NOT NULL,

  PRIMARY KEY (relation_id),
  UNIQUE INDEX child_parent_idx
    (child_term_id, parent_term_id, relation_type_id),
  INDEX parent_idx (parent_term_id)
);

DROP TABLE IF EXISTS closure;
CREATE TABLE closure (
  closure_id        INT UNSIGNED NOT NULL AUTO_INCREMENT,
  child_term_id     INT UNSIGNED NOT NULL,
  parent_term_id    INT UNSIGNED NOT NULL,
  subparent_term_id INT UNSIGNED,
  distance          TINYINT UNSIGNED NOT NULL,

  PRIMARY KEY (closure_id),
  UNIQUE INDEX child_parent_idx
    (child_term_id, parent_term_id, subparent_term_id),
  INDEX parent_subparent_idx
    (parent_term_id, subparent_term_id)
);

-- There are additional tables in the released databases called
-- "aux_XX_YY_map".  These are created by the "add_subset_maps.pl"
-- scripts.  Please see the README document for further information.

-- $Id$
