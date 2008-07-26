# patch_50_51_f.sql
#
# title: Set meta species_id
#
# description:
# Set the species_id of non-species-specific rows in the meta table to be NULL.

UPDATE  meta SET species_id = NULL WHERE meta_key IN ('patch', 'schema_version');

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_50_51_f.sql|meta_species_id_values');


