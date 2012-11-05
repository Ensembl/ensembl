# patch_69_70_d.sql
#
# Title: Restore data_file AUTO_INCREMENT field
#
# Description: patch_69_70_c.sql erased the AUTO_INCREMENT from data_file_id. 
#              This patch brings it back.

ALTER TABLE data_file MODIFY COLUMN data_file_id INT(10) UNSIGNED NOT NULL AUTO_INCREMENT;

# Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_69_70_d.sql|data_file_id_auto_increment');


