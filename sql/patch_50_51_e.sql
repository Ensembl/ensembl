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

# patch_50_51_e.sql
#
# title: Add external_data column to protein_feature and dna_align_feature.
#
# description:
# Add external_data column to protein_feature and dna_align_feature tables, primarily for storage of user-uploaded data.

ALTER TABLE protein_feature ADD COLUMN external_data TEXT;

ALTER TABLE dna_align_feature ADD COLUMN external_data TEXT;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_50_51_e.sql|feature_external_data');


