# patch_48_49_a.sql
#
# title: update schema version
#
# description:
# Update schema_version in meta table to 49

UPDATE meta SET meta_value='49' WHERE meta_key='schema_version';

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_48_49_a.sql|schema_version');


