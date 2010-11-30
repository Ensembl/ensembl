-- Schema of tables not added by the bootstrap_master.pl script.

-- NB: Additional tables are added by the web team to support storing
-- declarations of intentions etc.  Those tables are not even mentioned
-- here.

-- The 'species' table.
-- Lists the species for which there is a Core database.
CREATE TABLE species (
  species_id    INTEGER UNSIGNED NOT NULL AUTO_INCREMENT,
  db_name       VARCHAR(32) NOT NULL,   -- Name used in database names.
  common_name   VARCHAR(32) NOT NULL,   -- What we often refer to it as.
  web_name      VARCHAR(32) NOT NULL,   -- Name that the web site is using.
  is_current    BOOLEAN NOT NULL DEFAULT true,

  PRIMARY KEY (species_id),
  UNIQUE INDEX db_name_idx (db_name)
);


-- The 'db' table.
-- This table contains all species-specific databases for this release.
CREATE TABLE db (
  db_id         INTEGER UNSIGNED NOT NULL AUTO_INCREMENT,
  species_id    INTEGER UNSIGNED NOT NULL,  -- FK into 'species'.
  db_type       ENUM('cdna', 'core', 'coreexpressionatlas',
                     'coreexpressionest', 'coreexpressiongnf',
                     'funcgen', 'otherfeatures', 'variation', 'vega')
                     NOT NULL DEFAULT 'core',
  db_release    INTEGER NOT NULL,
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
                    'otherfeatures', 'variation', 'vega')
                    NOT NULL DEFAULT 'core',
  description   TEXT,

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
                      'variation', 'vega') NOT NULL DEFAULT 'core',
  only_for_species  TEXT,
  description       TEXT,

  PRIMARY KEY (meta_key_id),
  UNIQUE INDEX name_type_idx (name, db_type)
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

  PRIMARY KEY (analysis_description_id),
  UNIQUE INDEX logic_name_idx (logic_name)
);

-- The 'web_data' table.
-- TODO: ANY DATA FOUND IN THIS TABLE IS NOT YET "REAL".
--       DEVELOPMENT IS STILL UNDERWAY.
-- Contains the data for the 'web_data' column in the
-- 'analysis_description' table.
-- The 'web_data' is a hash and we store this as key-value pairs
-- ('hash_key' and 'hash_value').  The 'hash_key' might contain double
-- colons ('::') to distinguish sub-hash keys, e.g. 'default::MultiTop'.
CREATE TABLE web_data (
  web_data_id   INTEGER UNSIGNED NOT NULL AUTO_INCREMENT,
  hash_key      VARCHAR(32) NOT NULL,
  hash_value    VARCHAR(128),

  PRIMARY KEY (web_data_id)
);

-- The 'analysis_web_data' table.
-- TODO: ANY DATA FOUND IN THIS TABLE IS NOT YET "REAL".
--       DEVELOPMENT IS STILL UNDERWAY.
-- This table connects the 'analysis_description' table with the
-- 'web_data' and 'db' tables.
CREATE TABLE analysis_web_data (
  analysis_description_id   INTEGER UNSIGNED NOT NULL,
  web_data_id               INTEGER UNSIGNED NOT NULL,
  db_id                     INTEGER UNSIGNED NOT NULL,

  UNIQUE KEY analysis_web_data_db_idx
    (analysis_description_id, web_data_id, db_id)
);

-- VIEWS

CREATE VIEW db_list AS
SELECT  db_id AS db_id,
        CONCAT(
          CONCAT_WS('_', db_name, db_type, db_release, db_assembly),
        db_suffix) AS full_db_name
FROM    species
  JOIN  db USING (species_id);

-- CREATE VIEW readable_web_data AS
-- SELECT  CONCAT('{',
--           GROUP_CONCAT(data SEPARATOR ','),
--         '}') AS web_data
-- FROM    analysis_web_data awd
--   JOIN  web_data wd USING (web_data_id)
-- GROUP BY    
