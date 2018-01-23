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

# patch_69_70_b.sql
#
# Title: Update mapping_set to store history over several releases
#
# Description:

# Updating schema_build row to current_schema_build
# Add new row old_schema_build
# Remove row schema_build, to be replaced by current_schema_build
# All previously stored data will be deleted

ALTER TABLE mapping_set DROP COLUMN schema_build;
ALTER TABLE mapping_set ADD COLUMN (internal_schema_build VARCHAR(20) NOT NULL);
ALTER TABLE mapping_set ADD COLUMN (external_schema_build VARCHAR(20) NOT NULL);
TRUNCATE TABLE mapping_set;
ALTER TABLE mapping_set ADD UNIQUE KEY mapping_idx (internal_schema_build, external_schema_build);
TRUNCATE TABLE seq_region_mapping;


# Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_69_70_b.sql|add_mapping_set_history');


