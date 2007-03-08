# patch_43_44_e
#
# title: update schema version
#
# description:
# Update schema_version in meta table to 44

UPDATE meta SET meta_value='44' WHERE meta_key='schema_version';

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_43_44_e.sql|schema_version');


