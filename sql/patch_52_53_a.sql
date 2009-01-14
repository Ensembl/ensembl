# patch_52_53_a.sql
#
# title: update schema version
#
# description:
# Update schema_version in meta table to 53

UPDATE meta SET meta_value='53' WHERE meta_key='schema_version';

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_52_53_a.sql|schema_version');


