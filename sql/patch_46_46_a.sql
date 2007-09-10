# patch_46_47_a.sql
#
# title: update schema version
#
# description:
# Update schema_version in meta table to 47

UPDATE meta SET meta_value='47' WHERE meta_key='schema_version';

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_46_47_a.sql|schema_version');


