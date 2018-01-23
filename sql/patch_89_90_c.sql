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

# patch_89_90_c.sql
#
# Title: Protein_feature hit_name case sensitive.
#
# Description:
#   Make hit_name in protein_feature table case sensitive.
#   Required to store PDB features like 'abcd.a' and 'abcd.A'.

ALTER TABLE protein_feature CHANGE hit_name hit_name VARCHAR(40) CHARACTER SET latin1 COLLATE latin1_bin NOT NULL;

INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_89_90_c.sql|pf_hit_name_case_sensitive');
