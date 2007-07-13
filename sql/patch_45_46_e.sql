# patch_45_46_e.sql
#
# title: External_db new columns
#
# description:
# Add secondary_db_name/table columns to external_db - required for integration with eFG

ALTER TABLE external_db ADD COLUMN secondary_db_name  VARCHAR(255) DEFAULT NULL;
ALTER TABLE external_db ADD COLUMN secondary_db_table VARCHAR(255) DEFAULT NULL;

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_45_46_e.sql|external_db_new_cols');
