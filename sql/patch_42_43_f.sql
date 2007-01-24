# patch_42_43_f
#
# title: update schema version
#
# description:
# Update schema_version in meta table to 43

UPDATE meta SET meta_value='43' WHERE meta_key='schema_version';

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_42_43_f.sql|schema_version');


