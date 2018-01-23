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

# patch_63_64_c.sql
#
# Title:
#   is_ref to alt_allele 
#
# Description:
#    Add is_ref to alt_allele to specify which is the reference al


# Add the new column.
ALTER TABLE alt_allele
  ADD COLUMN is_ref BOOLEAN NOT NULL DEFAULT 0
    AFTER gene_id;

# Insert patch identifier.
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch',
    'patch_63_64_c.sql|is_ref_added_to_alt_allele');
