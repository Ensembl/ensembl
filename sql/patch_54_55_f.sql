# patch_54_55_f.sql
#
# title: Set to NULL version column in coord_system table instead of blank.
#
# description:
# The version column in the coord_system table has to be NULL when there
# is no version for the coordinate system (at the moment, some had blank entries)

UPDATE coord_system SET version = NULL WHERE version = '';

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
VALUES (NULL, 'patch', 'patch_54_55_f.sql|coord_system.version_null');
