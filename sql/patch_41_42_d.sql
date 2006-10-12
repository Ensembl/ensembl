# patch_40_41_d
#
# title: update schema version
#
# description:
# Update schema_version in meta table to 42

UPDATE meta SET meta_value='42' WHERE meta_key='schema_version';

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_41_42_d.sql|schema_version');


