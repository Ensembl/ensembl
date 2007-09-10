# patch_46_47_e.sql
#
# title: auto_increment external_db_id
#
# description:
# Add AUTO_INCREMENT to external_db.external_db_id

ALTER table external_db CHANGE COLUMN external_db_id external_db_id SMALLINT UNSIGNED NOT NULL AUTO_INCREMENT;

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_46_47_e.sql|auto_increment_external_db');


