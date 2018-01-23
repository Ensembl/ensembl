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

# patch_85_86_b.sql
#
# Title: New dna_align_feature_attrib table
#
# Description:
#   Add a new table, for attributes associated with DNA alignment features

CREATE TABLE dna_align_feature_attrib (

  dna_align_feature_id        INT(10) UNSIGNED NOT NULL,
  attrib_type_id              SMALLINT(5) UNSIGNED NOT NULL,
  value                       TEXT NOT NULL,

  UNIQUE KEY dna_align_feature_attribx (dna_align_feature_id, attrib_type_id, value(500)),
  KEY dna_align_feature_idx (dna_align_feature_id),
  KEY type_val_idx (attrib_type_id, value(40)),
  KEY val_only_idx (value(40))

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

# Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_85_86_b.sql|add dna_align_feature_attrib table');
