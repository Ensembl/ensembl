# patch_50_51_h.sql
#
# title: External_db name
#
# description:
# Extend db_name in external_db table to 40 characters (from 28)

ALTER TABLE external_db CHANGE COLUMN db_name db_name VARCHAR(40) NOT NULL;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_50_51_h.sql|external_db_db_name');


