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

# patch_83_84_e.sql
#
# Title: Make stable ID version nullable
#
# Description:
#   Make version nullable for gene, transcript, translation, exon, operon and operon_transcript

ALTER TABLE gene MODIFY COLUMN version SMALLINT UNSIGNED DEFAULT NULL;
ALTER TABLE transcript MODIFY COLUMN version SMALLINT UNSIGNED DEFAULT NULL;
ALTER TABLE translation MODIFY COLUMN version SMALLINT UNSIGNED DEFAULT NULL;
ALTER TABLE exon MODIFY COLUMN version SMALLINT UNSIGNED DEFAULT NULL;
ALTER TABLE operon MODIFY COLUMN version SMALLINT UNSIGNED DEFAULT NULL;
ALTER TABLE operon_transcript MODIFY COLUMN version SMALLINT UNSIGNED DEFAULT NULL;

# Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_83_84_e.sql|nullable_versions');
