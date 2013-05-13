-- patch_71_72c.sql
--
-- Title: Insert schema version.
--
-- Description:
--   Adding schema version to the meta table (set to 72)
--   so that script schema_patcher.pl would work

INSERT INTO meta (meta_key, meta_value)
  VALUES ('schema_version', 72);

-- Patch identifier
INSERT INTO meta (meta_key, meta_value)
  VALUES ('patch', 'patch_71_72c.sql|schema_version');


