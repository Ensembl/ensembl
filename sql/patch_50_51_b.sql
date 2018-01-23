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

# patch_50_51_b.sql
#
# title: protein_feature hit_name
#
# description:
# Rename hit_id column in protein_feature to hit_name

ALTER TABLE protein_feature DROP INDEX hid_index;

ALTER TABLE protein_feature CHANGE COLUMN hit_id hit_name VARCHAR(40) NOT NULL;

ALTER TABLE protein_feature ADD INDEX hitname_idx (hit_name);

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_50_51_b.sql|protein_feature_hit_name');


