-- patch_71_72_d.sql
--
-- Title: Fix patch versions
--
-- Description:
--   Fixes the existing patch meta items as their versioning was not correct

update meta set meta_value = 'patch_71_72_b.sql|alt_id table'
  where meta_key = 'patch'
  and meta_value = 'patch_71_72b.sql|alt_id table';

update meta set meta_value = 'patch_71_72_c.sql|schema_version'
  where meta_key = 'patch'
  and meta_value = 'patch_71_72c.sql|schema_version';

-- Patch identifier
INSERT INTO meta (meta_key, meta_value)
  VALUES ('patch', 'patch_71_72_d.sql|patch_version_fix');


