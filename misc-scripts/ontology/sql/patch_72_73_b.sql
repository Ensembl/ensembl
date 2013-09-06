-- patch_72_73_a.sql
--
-- Title: Insert schema version.
--
-- Description:
--   Adding schema version to the meta table (set to 73)
--   so that script schema_patcher.pl would work

ALTER TABLE meta ADD COLUMN species_id INT(1) UNSIGNED DEFAULT NULL;

-- Patch identifier
INSERT INTO meta (meta_key, meta_value)
  VALUES ('patch', 'patch_72_73_b.sql|meta_species');


