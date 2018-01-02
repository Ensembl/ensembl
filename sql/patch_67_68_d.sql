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

# patch_67_68_d.sql
#
# Title: Altering intron_supporting_evidence table
#
# Description: Adding is_splice_canonical which lets us define if a splice 
#              junction can be considered canonical and adding missing indexes
# 

ALTER TABLE intron_supporting_evidence 
ADD COLUMN is_splice_canonical BOOLEAN NOT NULL DEFAULT 0;

ALTER TABLE intron_supporting_evidence ADD KEY seq_region_idx (seq_region_id, seq_region_start);

# Patch identifier:
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_67_68_d.sql|add_is_splice_canonical_and_seq_index');
