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

# patch_56_57_f.sql
#
# Title: Expand width of simple_feature.display_label
#
# Description:
# The display_label field of the simple_feature table needs to be
# expanded from 40 characters to be able to hold some of the data
# imported from Encode for the human database.

ALTER TABLE simple_feature
MODIFY COLUMN display_label VARCHAR(255) NOT NULL;

OPTIMIZE TABLE simple_feature;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
VALUES (NULL, 'patch', 'patch_56_57_f.sql|simple_feature.display_label');
