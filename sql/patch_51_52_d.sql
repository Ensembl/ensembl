# patch_51_52_d.sql
#
# title: External_db description
#
# description:
# Add a description column to external_db

ALTER TABLE external_db ADD COLUMN description TEXT;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_51_52_d.sql|external_db_description');


