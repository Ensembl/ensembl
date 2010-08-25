# patch_59_60_b.sql
#
# Title:
#   Rename 'go_xref' table to 'ontology_xref'.
#
# Description:
#   Rename the 'go_xref' table to make its use more generic.

# Rename the table, and swap the source_xref_id and linkage_type fields.
ALTER TABLE go_xref
  RENAME TO ontology_xref,
  MODIFY COLUMN source_xref_id INT(10) UNSIGNED DEFAULT NULL
    AFTER object_xref_id;

# Optimize the table, because indexes may be out of whack
OPTIMIZE TABLE ontology_xref;

# Insert patch identifier.
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_59_60_b.sql|rename_go_xref_table');
