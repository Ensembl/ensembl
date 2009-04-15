# patch_51_52_a.sql
#
# title: update schema version
#
# description:
# Update schema_version in meta table to 52

UPDATE meta SET meta_value='52' WHERE meta_key='schema_version';

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_51_52_a.sql|schema_version');


