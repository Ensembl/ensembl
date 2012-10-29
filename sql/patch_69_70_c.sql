# patch_69_70_c.sql
#
# Title: Ensure column definitions are consistent in the schema
#
# Description: A number of column defintions over time have diverged from their
#              original specification. We are converting those we know are wrong

ALTER TABLE dependent_xref MODIFY COLUMN object_xref_id INT(10) UNSIGNED NOT NULL;
ALTER TABLE dependent_xref MODIFY COLUMN master_xref_id INT(10) UNSIGNED NOT NULL;
ALTER TABLE dependent_xref MODIFY COLUMN dependent_xref_id INT(10) UNSIGNED NOT NULL;

ALTER TABLE object_xref MODIFY COLUMN xref_id INT(10) UNSIGNED NOT NULL;

ALTER TABLE data_file MODIFY COLUMN data_file_id INT(10) UNSIGNED NOT NULL;
ALTER TABLE data_file MODIFY COLUMN coord_system_id INT(10) UNSIGNED NOT NULL;
ALTER TABLE data_file MODIFY COLUMN analysis_id SMALLINT UNSIGNED NOT NULL;

ALTER TABLE data_file COLLATE=latin1_swedish_ci;

# Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_69_70_c.sql|column_datatype_consistency');


