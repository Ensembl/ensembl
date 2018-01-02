-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2018] EMBL-European Bioinformatics Institute
-- 
-- Licensed under the Apache License, Version 2.0 (the "License");
-- you may not use this file except in compliance with the License.
-- You may obtain a copy of the License at
-- 
--      http://www.apache.org/licenses/LICENSE-2.0
-- 
-- Unless required by applicable law or agreed to in writing, software
-- distributed under the License is distributed on an "AS IS" BASIS,
-- WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
-- See the License for the specific language governing permissions and
-- limitations under the License.

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


