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

# patch_46_47_b.sql
#
# title: new align columns
#
# description:
# Add external_db_id and hcoverage columns to dna_align_feature and protein_align_feature

ALTER TABLE dna_align_feature ADD COLUMN external_db_id SMALLINT UNSIGNED;
ALTER TABLE dna_align_feature ADD COLUMN hcoverage DOUBLE;

ALTER TABLE dna_align_feature ADD KEY external_db_idx (external_db_id);

ALTER TABLE protein_align_feature ADD COLUMN external_db_id SMALLINT UNSIGNED;
ALTER TABLE protein_align_feature ADD COLUMN hcoverage DOUBLE;

ALTER TABLE protein_align_feature ADD KEY external_db_idx (external_db_id);

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_46_47_b.sql|new_align_columns');


