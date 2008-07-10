# patch_49_50_a.sql
#
# title: update schema version
#
# description:
# Update schema_version in meta table to 50

UPDATE meta SET meta_value='50' WHERE meta_key='schema_version';

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_49_50_a.sql|schema_version');


