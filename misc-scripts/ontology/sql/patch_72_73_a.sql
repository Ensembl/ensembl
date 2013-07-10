-- patch_72_73_a.sql
--
-- Title: Insert schema version.
--
-- Description:
--   Adding schema version to the meta table (set to 73)
--   so that script schema_patcher.pl would work

INSERT INTO meta (meta_key, meta_value)
  VALUES ('schema_version', 73);

-- Patch identifier
INSERT INTO meta (meta_key, meta_value)
  VALUES ('patch', 'patch_72_73_a.sql|schema_version');


