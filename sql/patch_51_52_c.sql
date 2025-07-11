-- See the NOTICE file distributed with this work for additional information
-- regarding copyright ownership.
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

# patch_51_52_c.sql
#
# title: Add dna_align_feature.pair_dna_align_feature_id
#
# description:
# Add a pair_dna_align_feature_id to allow support for paired reads

ALTER TABLE dna_align_feature ADD COLUMN pair_dna_align_feature_id INT(10) UNSIGNED;

ALTER TABLE dna_align_feature ADD INDEX pair_idx (pair_dna_align_feature_id);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_51_52_c.sql|pair_dna_align_feature_id');


