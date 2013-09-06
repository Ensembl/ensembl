# patch_73_74_a.sql
#
# Title: Update schema version.
#
# Description:
#   Update schema_version in meta table to 74.

UPDATE meta SET meta_value='74' WHERE meta_key='schema_version';

# Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_73_74_a.sql|schema_version');
