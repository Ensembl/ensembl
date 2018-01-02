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

# patch_60_61_b.sql
#
# Title:
#   Create 'seq_region_synonym' table.
#
# Description: Create 'seq_region_synonym' table which will hold the alternative seq region 
# names.

# Creat the tab
CREATE TABLE seq_region_synonym (

  seq_region_synonym_id       INT UNSIGNED NOT NULL  AUTO_INCREMENT,
  seq_region_id               INT(10) UNSIGNED NOT NULL,
  synonym                     VARCHAR(40) NOT NULL,
  external_db_id              SMALLINT UNSIGNED,

  PRIMARY KEY (seq_region_synonym_id),
  UNIQUE KEY syn_idx (synonym)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

# Insert patch identifier.
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_60_61_b.sql|create_seq_region_synonym_table');
