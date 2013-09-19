-- patch_73_74_a.sql
--
-- Title: Insert schema version.
--
-- Description:
--   Adding schema version to the meta table (set to 74)
--   so that script schema_patcher.pl would work

INSERT INTO meta (meta_key, meta_value)
  VALUES ('schema_version', 74);

-- Patch identifier
INSERT INTO meta (meta_key, meta_value)
  VALUES ('patch', 'patch_73_74_a.sql|schema_version');


