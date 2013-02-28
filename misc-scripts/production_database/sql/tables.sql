-- Table schema for the production database


-- The 'species' table.
-- Lists the species for which there is a Core database.
CREATE TABLE species (
  species_id      INTEGER UNSIGNED NOT NULL AUTO_INCREMENT,
  db_name         VARCHAR(255) NOT NULL, -- Name used in database names.
  common_name     VARCHAR(255) NOT NULL, -- What we often refer to it as.
  web_name        VARCHAR(255) NOT NULL, -- Name that the web site is using.
  scientific_name VARCHAR(255) NOT NULL, -- Full name of the species
  production_name VARCHAR(255) NOT NULL, -- Name that production processes use
  url_name        VARCHAR(255) NOT NULL, -- Name that is used in URLs
  taxon           VARCHAR(20),
  species_prefix  VARCHAR(20),
  is_current      BOOLEAN NOT NULL DEFAULT true,
  attrib_type_id  SMALLINT(5) UNSIGNED DEFAULT NULL,

  -- Columns for the web interface:
  created_by    INTEGER,
  created_at    DATETIME,
  modified_by   INTEGER,
  modified_at   DATETIME,

  PRIMARY KEY (species_id),
  UNIQUE INDEX db_name_idx (db_name)
);

-- The 'species_alias' table
-- Lists all aliases for all species including those no longer active

CREATE TABLE species_alias (
  species_alias_id  INTEGER UNSIGNED NOT NULL AUTO_INCREMENT, -- surrogate key
  species_id        INTEGER UNSIGNED NOT NULL,      -- FK into species
  alias             varchar(255) NOT NULL,          -- alias
  is_current        BOOLEAN NOT NULL DEFAULT true,  -- if it's still current
  
  -- Columns for the web interface:
  created_by    INTEGER,
  created_at    DATETIME,
  modified_by   INTEGER,
  modified_at   DATETIME,
  
  PRIMARY KEY (species_alias_id),
  UNIQUE INDEX (alias, is_current),             -- aliases MUST be unique for 
                                                -- the current set. A certain 
                                                -- amount of duplication is 
                                                -- allowed if an alias moved 
                                                -- once
  INDEX sa_speciesid_idx (species_id)
);


-- The 'db' table.
-- This table contains all species-specific databases for this release.
CREATE TABLE db (
  db_id         INTEGER UNSIGNED NOT NULL AUTO_INCREMENT,
  species_id    INTEGER UNSIGNED NOT NULL,  -- FK into 'species'.
  is_current    BOOLEAN NOT NULL DEFAULT false,
  db_type       ENUM('cdna', 'core', 'coreexpressionatlas',
                     'coreexpressionest', 'coreexpressiongnf',
                     'funcgen', 'otherfeatures', 'rnaseq',
                     'variation', 'vega')
                     NOT NULL DEFAULT 'core',
  db_release    VARCHAR(8) NOT NULL,
  db_assembly   VARCHAR(8) NOT NULL,
  db_suffix     CHAR(1) DEFAULT '',
  db_host       VARCHAR(32) DEFAULT NULL,

  PRIMARY KEY (db_id),
  UNIQUE INDEX species_release_idx (species_id, db_type, db_release)
);


-- The 'biotype' table.
-- Contains all the valid biotypes used for genes and transcripts.
CREATE TABLE biotype (
  biotype_id    INTEGER UNSIGNED NOT NULL AUTO_INCREMENT,
  name          VARCHAR(64) NOT NULL,
  is_current    BOOLEAN NOT NULL DEFAULT true,
  is_dumped     BOOLEAN NOT NULL DEFAULT true,
  object_type   ENUM('gene', 'transcript') NOT NULL DEFAULT 'gene',
  db_type       SET('cdna', 'core', 'coreexpressionatlas',
                    'coreexpressionest', 'coreexpressiongnf', 'funcgen',
                    'otherfeatures', 'rnaseq', 'variation', 'vega')
                    NOT NULL DEFAULT 'core',
  description   TEXT,

  -- Columns for the web interface:
  created_by    INTEGER,
  created_at    DATETIME,
  modified_by   INTEGER,
  modified_at   DATETIME,

  PRIMARY KEY (biotype_id),
  UNIQUE INDEX name_type_idx (name, object_type, db_type)
);

-- The 'meta_key' table.
-- Contains the meta keys that may or must be available in the 'meta'
-- table in the Core databases.
CREATE TABLE meta_key (
  meta_key_id       INTEGER UNSIGNED NOT NULL AUTO_INCREMENT,
  name              VARCHAR(64) NOT NULL,
  is_optional       BOOLEAN NOT NULL DEFAULT false,
  is_current        BOOLEAN NOT NULL DEFAULT true,
  db_type           SET('cdna', 'core', 'funcgen', 'otherfeatures',
                        'rnaseq', 'variation', 'vega')
                    NOT NULL DEFAULT 'core',
  description       TEXT,

  -- Columns for the web interface:
  created_by    INTEGER,
  created_at    DATETIME,
  modified_by   INTEGER,
  modified_at   DATETIME,

  PRIMARY KEY (meta_key_id),
  UNIQUE INDEX name_type_idx (name, db_type)
);

-- The 'meta_key_species' table.
-- Connects the 'meta_key' and the 'species' tables.
CREATE TABLE meta_key_species (
  meta_key_id       INTEGER UNSIGNED NOT NULL,
  species_id        INTEGER UNSIGNED NOT NULL,

  UNIQUE INDEX uniq_idx (meta_key_id, species_id)
);

-- The 'analysis_description' table.
-- Contains the analysis logic name along with the data that should
-- be available in the 'analysis_description' table, except for the
-- 'web_data' and 'displayable' columns.
CREATE TABLE analysis_description (
  analysis_description_id   INTEGER UNSIGNED NOT NULL AUTO_INCREMENT,
  logic_name                VARCHAR(128) NOT NULL,
  description               TEXT,
  display_label             VARCHAR(256) NOT NULL,
  db_version                TINYINT(1) NOT NULL DEFAULT '1',
  is_current                BOOLEAN NOT NULL DEFAULT true,
  default_web_data_id       INT(10) UNSIGNED DEFAULT NULL,

  -- Columns for the web interface:
  created_by    INTEGER,
  created_at    DATETIME,
  modified_by   INTEGER,
  modified_at   DATETIME,

  PRIMARY KEY (analysis_description_id),
  UNIQUE INDEX logic_name_idx (logic_name)
);

-- The 'analysis_web_data' table.
-- Many-to-many connection table.
-- Contains the data for the 'displayable' columns in the
-- 'analysis_description' table.  Ties together species,
-- analysis_description, and the web_data.
CREATE TABLE analysis_web_data (
  analysis_web_data_id      INTEGER UNSIGNED NOT NULL AUTO_INCREMENT,
  analysis_description_id   INTEGER UNSIGNED NOT NULL,
  web_data_id               INTEGER UNSIGNED DEFAULT NULL,
  species_id                INTEGER UNSIGNED NOT NULL,

  db_type                   ENUM('cdna', 'core', 'funcgen',
                                'otherfeatures', 'rnaseq', 'vega',
                                'presite', 'sangervega')
                            NOT NULL DEFAULT 'core',

  displayable               BOOLEAN NOT NULL DEFAULT true,

  -- Columns for the web interface:
  created_by    INTEGER,
  created_at    DATETIME,
  modified_by   INTEGER,
  modified_at   DATETIME,

  PRIMARY KEY (analysis_web_data_id),
  UNIQUE INDEX uniq_idx (species_id, db_type, analysis_description_id),
  INDEX ad_idx (analysis_description_id),
  INDEX wd_idx (web_data_id)
);

-- The 'web_data' table.
-- Contains the unique web_data.
CREATE TABLE web_data (
  web_data_id               INTEGER UNSIGNED NOT NULL AUTO_INCREMENT,
  data                      TEXT,
  
  -- Columns for internal documentation
  comment TEXT,

  -- Columns for the web interface:
  created_by    INTEGER,
  created_at    DATETIME,
  modified_by   INTEGER,
  modified_at   DATETIME,

  PRIMARY KEY (web_data_id)
);

-- The 'master_attrib_type' table.
-- Lists the existing attrib_types
CREATE TABLE master_attrib_type (
  attrib_type_id            SMALLINT(5) unsigned NOT NULL AUTO_INCREMENT,
  code                      VARCHAR(15) NOT NULL DEFAULT '',
  name                      VARCHAR(255) NOT NULL DEFAULT '',
  description               TEXT,
  is_current                TINYINT(1) NOT NULL DEFAULT '1',

  -- Columns for the web interface:
  created_by    INTEGER,
  created_at    DATETIME,
  modified_by   INTEGER,
  modified_at   DATETIME,

  PRIMARY KEY (attrib_type_id),
  UNIQUE INDEX code_idx (code)
);

-- The 'master_external_db' table.
-- Lists the existing external_dbs
CREATE TABLE master_external_db (
  external_db_id            INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  db_name                   VARCHAR(100) NOT NULL,
  db_release                VARCHAR(255) DEFAULT NULL,
  status                    ENUM('KNOWNXREF','KNOWN','XREF','PRED','ORTH','PSEUDO') NOT NULL,
  priority                  INT(11) NOT NULL,
  db_display_name           VARCHAR(255) DEFAULT NULL,
  type                      ENUM('ARRAY','ALT_TRANS','ALT_GENE','MISC','LIT','PRIMARY_DB_SYNONYM','ENSEMBL') DEFAULT NULL,
  secondary_db_name         VARCHAR(255) DEFAULT NULL,
  secondary_db_table        VARCHAR(255) DEFAULT NULL,
  description               TEXT,
  is_current                TINYINT(1) NOT NULL DEFAULT '1',

  -- Columns for the web interface:
  created_by    INTEGER,
  created_at    DATETIME,
  modified_by   INTEGER,
  modified_at   DATETIME,

  PRIMARY KEY (external_db_id),
  UNIQUE INDEX db_name_idx (db_name,db_release,is_current)
);

-- The 'master_misc_set' table
-- Lists the existing misc_sets
CREATE TABLE master_misc_set (
  misc_set_id               SMALLINT(5) UNSIGNED NOT NULL AUTO_INCREMENT,
  code                      VARCHAR(25) NOT NULL DEFAULT '',
  name                      VARCHAR(255) NOT NULL DEFAULT '',
  description               TEXT NOT NULL,
  max_length                INT(10) UNSIGNED NOT NULL,
  is_current                TINYINT(1) NOT NULL DEFAULT '1',

  -- Columns for the web interface:
  created_by    INTEGER,
  created_at    DATETIME,
  modified_by   INTEGER,
  modified_at   DATETIME,

  PRIMARY KEY (misc_set_id),
  UNIQUE INDEX code_idx (code)
);


-- The 'master_unmapped_reason' table
-- Lists the existing unmapped_reasons
CREATE TABLE master_unmapped_reason (
  unmapped_reason_id        SMALLINT(5) UNSIGNED NOT NULL AUTO_INCREMENT,
  summary_description       VARCHAR(255) DEFAULT NULL,
  full_description          VARCHAR(255) DEFAULT NULL,
  is_current                TINYINT(1) NOT NULL DEFAULT '1',

  -- Columns for the web interface:
  created_by    INTEGER,
  created_at    DATETIME,
  modified_by   INTEGER,
  modified_at   DATETIME,

  PRIMARY KEY (unmapped_reason_id)
);


-- The 'changelog' and 'changelog_species' tables
-- Used for web purposes to store release related information
CREATE TABLE changelog (
  changelog_id              INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  release_id                INT(11) DEFAULT NULL,
  title                     VARCHAR(128) DEFAULT NULL,
  content                   TEXT,
  notes                     TEXT,
  status                    ENUM('declared','handed_over','postponed','cancelled') NOT NULL DEFAULT 'declared',
  team                      ENUM('Compara','Core','Funcgen','Genebuild','Outreach','Variation','Web','EnsemblGenomes','Wormbase','Production') DEFAULT NULL,
  assembly                  ENUM('N','Y') NOT NULL DEFAULT 'N',
  gene_set                  ENUM('N','Y') NOT NULL DEFAULT 'N',
  repeat_masking            ENUM('N','Y') NOT NULL DEFAULT 'N',
  stable_id_mapping         ENUM('N','Y') NOT NULL DEFAULT 'N',
  affy_mapping              ENUM('N','Y') NOT NULL DEFAULT 'N',
  biomart_affected          ENUM('N','Y') NOT NULL DEFAULT 'N',
  variation_pos_changed     ENUM('N','Y') NOT NULL DEFAULT 'N',
  db_status                 ENUM('N/A','unchanged','patched','new') NOT NULL DEFAULT 'N/A',
  db_type_affected          SET('cdna','core','funcgen','otherfeatures','rnaseq','variation','vega') DEFAULT NULL,
  priority                  TINYINT(1) NOT NULL DEFAULT '2',
  is_current                TINYINT(1) NOT NULL DEFAULT '1',

  -- Columns for the web interface:
  created_by    INTEGER,
  created_at    DATETIME,
  modified_by   INTEGER,
  modified_at   DATETIME,

  PRIMARY KEY (changelog_id)
);

-- Connects the 'changelog' and the 'species' tables.
CREATE TABLE changelog_species (
  changelog_id              INT(11) NOT NULL DEFAULT '0',
  species_id                INT(11) NOT NULL DEFAULT '0',

  PRIMARY KEY (changelog_id,species_id)
);

-- VIEWS

-- The 'db_list' view provides the full database names for all databases
-- in the 'db' table that are current.
CREATE DEFINER = CURRENT_USER SQL SECURITY INVOKER VIEW db_list AS
SELECT  db_id AS db_id,
        CONCAT(
          CONCAT_WS('_', db_name, db_type, db_release, db_assembly),
        db_suffix) AS full_db_name
FROM    species
  JOIN  db USING (species_id)
WHERE species.is_current = true;

-- The 'full_analysis_description' view provids, for each database,
-- /nearly/ exactly what should go into the 'analysis_description'
-- table, apart from the fact that it uses 'analysis.logic_name' rather
-- than 'analysis.analysis_id'.
CREATE DEFINER = CURRENT_USER SQL SECURITY INVOKER VIEW full_analysis_description AS
SELECT  list.full_db_name AS full_db_name,
        ad.logic_name AS logic_name,
        ad.description AS description,
        ad.display_label AS display_label,
        awd.displayable AS displayable,
        wd.data AS web_data
FROM db_list list
  JOIN db USING (db_id)
  JOIN analysis_web_data awd
    ON ( db.species_id = awd.species_id
    AND  db.db_type = awd.db_type )
  JOIN analysis_description ad USING (analysis_description_id)
  LEFT JOIN web_data wd USING (web_data_id)
WHERE   db.is_current = true
  AND   ad.is_current = true;

-- The 'logic_name_overview' is a helper view for people trying to
-- make sense of the 'analysis_description', 'analysis_web_data', and
-- 'web_data' tables.
CREATE DEFINER = CURRENT_USER SQL SECURITY INVOKER VIEW logic_name_overview AS
SELECT
  awd.analysis_web_data_id AS analysis_web_data_id,
  ad.logic_name AS logic_name,
  ad.analysis_description_id AS analysis_description_id,
  s.db_name AS species,
  s.species_id AS species_id,
  awd.db_type AS db_type,
  wd.web_data_id AS web_data_id,
  awd.displayable AS displayable
FROM   analysis_description ad
  JOIN analysis_web_data awd USING (analysis_description_id)
  JOIN species s USING (species_id)
  LEFT JOIN web_data wd USING (web_data_id)
WHERE   s.is_current = true
  AND   ad.is_current = true;

-- The 'unconnected_analyses' view gives back the analyses from
-- 'analysis_description' that are unused in 'analysis_web_data'.
CREATE DEFINER = CURRENT_USER SQL SECURITY INVOKER VIEW unconnected_analyses AS
SELECT  ad.analysis_description_id AS analysis_description_id,
        ad.logic_name AS logic_name
FROM    analysis_description ad
  LEFT JOIN analysis_web_data awd USING (analysis_description_id)
WHERE   awd.species_id IS NULL
  AND   ad.is_current = true;

-- The 'unused_web_data' view gives back the entries in 'web_data' that
-- are not connected to any analysis in 'analysis_web_data'.
CREATE DEFINER = CURRENT_USER SQL SECURITY INVOKER VIEW unused_web_data AS
SELECT  wd.web_data_id
FROM    web_data wd
  LEFT JOIN analysis_web_data awd USING (web_data_id)
WHERE   awd.analysis_web_data_id IS NULL;


-- Views for the master tables.  These four views are simply selecting
-- the entries from the corresponding master table that have is_current
-- set to true.

CREATE DEFINER = CURRENT_USER SQL SECURITY INVOKER VIEW attrib_type AS
SELECT
  attrib_type_id AS attrib_type_id,
  code AS code,
  name AS name,
  description AS description
FROM    master_attrib_type
WHERE   is_current = true
ORDER BY attrib_type_id;

CREATE DEFINER = CURRENT_USER SQL SECURITY INVOKER VIEW external_db AS
SELECT
  external_db_id AS external_db_id,
  db_name AS db_name,
  db_release AS db_release,
  status AS status,
  priority AS priority,
  db_display_name AS db_display_name,
  type AS type,
  secondary_db_name AS secondary_db_name,
  secondary_db_table AS secondary_db_table,
  description AS description
FROM    master_external_db
WHERE   is_current = true
ORDER BY external_db_id;

CREATE DEFINER = CURRENT_USER SQL SECURITY INVOKER VIEW misc_set AS
SELECT
  misc_set_id AS misc_set_id,
  code AS code,
  name AS name,
  description AS description,
  max_length AS max_length
FROM    master_misc_set
WHERE   is_current = true
ORDER BY misc_set_id;

CREATE DEFINER = CURRENT_USER SQL SECURITY INVOKER VIEW unmapped_reason AS
SELECT
  unmapped_reason_id AS unmapped_reason_id,
  summary_description AS summary_description,
  full_description AS full_description
FROM    master_unmapped_reason
WHERE   is_current = true
ORDER BY unmapped_reason_id;
