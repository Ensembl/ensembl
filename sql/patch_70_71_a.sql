# patch_70_71_a.sql
#
# Title: Update schema version.
#
# Description:
#   Update schema_version in meta table to 71.

UPDATE meta SET meta_value='71' WHERE meta_key='schema_version';

# Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_70_71_a.sql|schema_version');
