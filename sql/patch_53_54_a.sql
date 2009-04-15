# patch_53_54_a.sql
#
# title: update schema version
#
# description:
# Update schema_version in meta table to 54

UPDATE meta SET meta_value='54' WHERE meta_key='schema_version';

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_53_54_a.sql|schema_version');


