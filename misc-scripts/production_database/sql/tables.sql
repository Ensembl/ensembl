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
  is_current                BOOLEAN NOT NULL DEFAULT true,

  PRIMARY KEY (analysis_description_id),
  UNIQUE INDEX logic_name_idx (logic_name)
);

-- The 'analysis_web_data' table.
-- Many-to-many connection table.
-- Contains the data for the 'displayable' columns in the
-- 'analysis_description' table.  Ties together species,
-- analysis_description, and the web_data.
CREATE TABLE analysis_web_data (
  analysis_description_id   INTEGER UNSIGNED NOT NULL,
  web_data_id               INTEGER UNSIGNED DEFAULT NULL,
  species_id                INTEGER UNSIGNED NOT NULL,

  db_type                   SET('cdna', 'core', 'funcgen',
                                'otherfeatures', 'rnaseq', 'vega')
                            NOT NULL DEFAULT 'core',

  displayable               BOOLEAN NOT NULL DEFAULT true,

  UNIQUE INDEX uniq_idx (species_id, db_type, analysis_description_id)
);

-- The 'web_data' table.
-- Contains the unique web_data.
CREATE TABLE web_data (
  web_data_id               INTEGER UNSIGNED NOT NULL AUTO_INCREMENT,
  data                      TEXT,

  PRIMARY KEY (web_data_id)
);

-- VIEWS

CREATE VIEW db_list AS
SELECT  db_id AS db_id,
        CONCAT(
          CONCAT_WS('_', db_name, db_type, db_release, db_assembly),
        db_suffix) AS full_db_name
FROM    species
  JOIN  db USING (species_id)
WHERE species.is_current = true;

CREATE VIEW full_analysis_description AS
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

CREATE VIEW logic_name_overview AS
SELECT
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

CREATE VIEW unconnected_analyses AS
SELECT  ad.analysis_description_id AS analysis_description_id,
        ad.logic_name AS logic_name
FROM    analysis_description ad
  LEFT JOIN analysis_web_data awd USING (analysis_description_id)
WHERE   awd.species_id IS NULL
  AND   ad.is_current = true;


-- Views for the master tables.  These are simply selecting the entries
-- from the corresponding master table that have is_current = 1.

CREATE VIEW attrib_type AS
SELECT
  attrib_type_id AS attrib_type_id,
  code AS code,
  name AS name,
  description AS description
FROM    master_attrib_type
WHERE   is_current = true;

CREATE VIEW external_db AS
SELECT
  external_db_id AS external_db_id,
  db_name AS db_name,
  db_release AS db_release,
  status AS status,
  dbprimary_acc_linkable AS dbprimary_acc_linkable,
  priority AS priority,
  db_display_name AS db_display_name,
  type AS type,
  secondary_db_name AS secondary_db_name,
  secondary_db_table AS secondary_db_table,
  description AS description
FROM    master_external_db
WHERE   is_current = true;

CREATE VIEW misc_set AS
SELECT
  misc_set_id AS misc_set_id,
  code AS code,
  name AS name,
  description AS description,
  max_length AS max_length
FROM    master_misc_set
WHERE   is_current = true;

CREATE VIEW unmapped_reason AS
SELECT
  unmapped_reason_id AS unmapped_reason_id,
  summary_description AS summary_description,
  full_description AS full_description
FROM    master_unmapped_reason
WHERE   is_current = true;
