# patch_55_56_a.sql
#
# title: update schema version
#
# description:
# Update schema_version in meta table to 56

UPDATE meta SET meta_value='56' WHERE meta_key='schema_version';

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_55_56_a.sql|schema_version');


