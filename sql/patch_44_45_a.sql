# patch_44_45_a.sql
#
# title: update schema version
#
# description:
# Update schema_version in meta table to 45

UPDATE meta SET meta_value='45' WHERE meta_key='schema_version';

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_44_45_a.sql|schema_version');


