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

# patch_86_87_c.sql
#
# Title: DATETIME DEFAULT NULL 
#
# Description:
# Later MySQLs no longer accept 0 in datetime fields and CURRENT_TIMESTAMP was introduced for dynamic default value.
# Older MySQLs do not support CURRENT_TIMESTAMP. To be able to support both, DATETIME DEFAULT is set to NULL

ALTER TABLE genome_statistics MODIFY COLUMN timestamp DATETIME NULL DEFAULT NULL;
ALTER TABLE analysis MODIFY COLUMN created DATETIME NULL DEFAULT NULL;
ALTER TABLE exon MODIFY COLUMN created_date DATETIME NULL DEFAULT NULL;
ALTER TABLE exon MODIFY COLUMN modified_date DATETIME NULL DEFAULT NULL;
ALTER TABLE gene MODIFY COLUMN created_date DATETIME NULL DEFAULT NULL;
ALTER TABLE gene MODIFY COLUMN modified_date DATETIME NULL DEFAULT NULL;
ALTER TABLE transcript MODIFY COLUMN created_date DATETIME NULL DEFAULT NULL;
ALTER TABLE transcript MODIFY COLUMN modified_date DATETIME NULL DEFAULT NULL;
ALTER TABLE translation MODIFY COLUMN created_date DATETIME NULL DEFAULT NULL;
ALTER TABLE translation MODIFY COLUMN modified_date DATETIME NULL DEFAULT NULL;
ALTER TABLE operon MODIFY COLUMN created_date DATETIME NULL DEFAULT NULL;
ALTER TABLE operon MODIFY COLUMN modified_date DATETIME NULL DEFAULT NULL;
ALTER TABLE operon_transcript MODIFY COLUMN created_date DATETIME NULL DEFAULT NULL;
ALTER TABLE operon_transcript MODIFY COLUMN modified_date DATETIME NULL DEFAULT NULL;

# Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_86_87_c.sql|datetime_default_NULL');
