# patch_68_69_a.sql
#
# Title: Update schema version.
#
# Description:
#   Update schema_version in meta table to 69.

UPDATE meta SET meta_value='69' WHERE meta_key='schema_version';

# Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_68_69_a.sql|schema_version');
