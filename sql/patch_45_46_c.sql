# patch_45_46_c.sql
#
# title: unmapped_object.external_db_id
#
# description:
# Make external_db column in unmapped_object be SMALLINT UNSIGNED explicitly.

ALTER TABLE unmapped_object CHANGE COLUMN external_db_id external_db_id SMALLINT UNSIGNED;

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_45_46_c.sql|unmapped_object.external_db_id');



