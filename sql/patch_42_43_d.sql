# patch_42_43_d
#
# title: Unmapped_object external_db_id
#
# description:
# Remove NOT NULL constraint on unmapped_object.external_db_id


ALTER TABLE unmapped_object CHANGE COLUMN external_db_id  external_db_id INT(10) UNSIGNED;

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_42_43_d.sql|unmapped_object_external_db_id');