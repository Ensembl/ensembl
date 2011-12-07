-- MAKING ALL VIEWS HAVE INVOKER SECURITY

CREATE OR REPLACE  DEFINER = CURRENT_USER SQL SECURITY INVOKER VIEW db_list AS
SELECT  db_id AS db_id,
        CONCAT(
          CONCAT_WS('_', db_name, db_type, db_release, db_assembly),
        db_suffix) AS full_db_name
FROM    species
  JOIN  db USING (species_id)
WHERE species.is_current = true;

CREATE OR REPLACE  DEFINER = CURRENT_USER SQL SECURITY INVOKER VIEW full_analysis_description AS
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

CREATE OR REPLACE  DEFINER = CURRENT_USER SQL SECURITY INVOKER VIEW logic_name_overview AS
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

CREATE OR REPLACE  DEFINER = CURRENT_USER SQL SECURITY INVOKER VIEW unconnected_analyses AS
SELECT  ad.analysis_description_id AS analysis_description_id,
        ad.logic_name AS logic_name
FROM    analysis_description ad
  LEFT JOIN analysis_web_data awd USING (analysis_description_id)
WHERE   awd.species_id IS NULL
  AND   ad.is_current = true;

CREATE OR REPLACE  DEFINER = CURRENT_USER SQL SECURITY INVOKER VIEW unused_web_data AS
SELECT  wd.web_data_id
FROM    web_data wd
  LEFT JOIN analysis_web_data awd USING (web_data_id)
WHERE   awd.analysis_web_data_id IS NULL;

CREATE OR REPLACE  DEFINER = CURRENT_USER SQL SECURITY INVOKER VIEW attrib_type AS
SELECT
  attrib_type_id AS attrib_type_id,
  code AS code,
  name AS name,
  description AS description
FROM    master_attrib_type
WHERE   is_current = true
ORDER BY attrib_type_id;

CREATE OR REPLACE  DEFINER = CURRENT_USER SQL SECURITY INVOKER VIEW external_db AS
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

CREATE OR REPLACE  DEFINER = CURRENT_USER SQL SECURITY INVOKER VIEW misc_set AS
SELECT
  misc_set_id AS misc_set_id,
  code AS code,
  name AS name,
  description AS description,
  max_length AS max_length
FROM    master_misc_set
WHERE   is_current = true
ORDER BY misc_set_id;

CREATE OR REPLACE  DEFINER = CURRENT_USER SQL SECURITY INVOKER VIEW unmapped_reason AS
SELECT
  unmapped_reason_id AS unmapped_reason_id,
  summary_description AS summary_description,
  full_description AS full_description
FROM    master_unmapped_reason
WHERE   is_current = true
ORDER BY unmapped_reason_id;
