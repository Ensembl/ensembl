# patch_56_57_d.sql
#
# Title: Allow for NULLs in meta table.
#
# Description:
# Allow the meta_value field in the meta table to be NULL.

ALTER TABLE meta
MODIFY COLUMN meta_value VARCHAR(255) BINARY;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
VALUES (NULL, 'patch', 'patch_56_57_d.sql|allow_meta_null');
