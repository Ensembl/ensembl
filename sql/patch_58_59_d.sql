# patch_58_59_d.sql
#
# Title:
#   Extend the object_type_idx index of the object_xref table.
#
# Description:
#   Add the 'analysis_id' to the end of the object_type_idx index in the
#   object_xref table so that we may allow for object_xref entries that
#   only differ in analysis_id.

# Modify the object_xref table.
ALTER TABLE object_xref
  DROP INDEX object_type_idx,
  ADD UNIQUE INDEX
    object_type_idx (ensembl_object_type, ensembl_id, xref_id, analysis_id);

ANALYZE TABLE object_xref;

# Insert patch identifier.
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch',
    'patch_58_59_d.sql|object_xref_extend_index');
