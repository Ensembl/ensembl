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


