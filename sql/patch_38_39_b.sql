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

# patch_38_39_b
#
# title: unique assembly
#
# description:
# this patch adds a unique index to the assembly table to prevent duplicate data

# create a new table with the unique index

CREATE TABLE assembly_new (

  asm_seq_region_id           INT UNSIGNED NOT NULL,
  cmp_seq_region_id           INT(10) UNSIGNED NOT NULL, 
  asm_start                   INT(10) NOT NULL,
  asm_end                     INT(10) NOT NULL,
  cmp_start                   INT(10) NOT NULL,
  cmp_end                     INT(10) NOT NULL,
  ori                         TINYINT  NOT NULL, 
  
  KEY (cmp_seq_region_id),
  KEY (asm_seq_region_id, asm_start),
  UNIQUE KEY all_idx (asm_seq_region_id, cmp_seq_region_id, asm_start, asm_end, cmp_start, cmp_end, ori)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;

# insert unique values into the new assembly table
INSERT IGNORE INTO assembly_new SELECT * FROM assembly;

# drop old assembly table and rename new one
DROP TABLE assembly;
ALTER TABLE assembly_new RENAME assembly;

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_38_39_b.sql|unique_assembly');

