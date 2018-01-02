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

# patch_66_67_c.sql
#
# Title: Adding intron_supporting_evidence table
#
# Description: Introns can be supported by an external feature. This gives a 
# weight to how much we believe the intron
# 

CREATE TABLE intron_supporting_evidence (
  intron_supporting_evidence_id INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  previous_exon_id INT(10) UNSIGNED NOT NULL,
  next_exon_id INT(10) UNSIGNED NOT NULL,
  hit_name VARCHAR(100) NOT NULL,
  score DECIMAL(10,3),
  score_type ENUM('NONE', 'DEPTH') DEFAULT 'NONE',
  
  PRIMARY KEY (intron_supporting_evidence_id),
  
  UNIQUE KEY (previous_exon_id, next_exon_id)
) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

# Patch identifier:
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_66_67_c.sql|adding_intron_supporting_evidence');
