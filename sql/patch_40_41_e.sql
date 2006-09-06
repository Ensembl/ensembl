# patch_40_41_e
#
# title: update schema version
#
# description:
# Update schema_version in meta table to 41

UPDATE meta SET meta_value='41' WHERE meta_key='schema_version';

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_40_41_e.sql|schema_version');


