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

# patch_81_82_c.sql
#
# Title: Update key for seq_region_synonym table
#
# Description:
#   Key for seq_region_synonym updated to include seq_region_id
#   For species with multiple assemblies, we can have the same synonym
#   for the same chromosome in different assemblies

ALTER TABLE seq_region_synonym 
   DROP INDEX syn_idx,
   ADD UNIQUE KEY syn_idx (synonym, seq_region_id);

# Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_81_82_c.sql|seq_synonym_key');
