# patch_58_59_e.sql
#
# Title:
#   Add new 'schema_type' meta key.
#
# Description:
#   All Ensembl databases with a 'meta' table should have a
#   'schema_type' key in them.  For Core datbases, the corresponding
#   value should be 'core'.

# Insert the new meta key.
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'schema_type', 'core');

# Insert patch identifier.
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_58_59_e.sql|meta_schema_type');
