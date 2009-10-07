# patch_56_57_a.sql
#
# title: update schema version
#
# description:
# Update schema_version in meta table to 57

UPDATE meta SET meta_value='57' WHERE meta_key='schema_version';

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_56_57_a.sql|schema_version');


