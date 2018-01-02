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

# patch_79_80_d.sql
#
# Title: Update column definition for genome_statistics.value from INT to BIGINT, 
#
# Description:
#   Column type from INT(10) to BIGINT(11) for genome_statistics.value
#   So value range can be extended to large plants genomes

ALTER TABLE genome_statistics MODIFY COLUMN value BIGINT(11) unsigned NOT NULL DEFAULT '0'; 

# Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_79_80_d.sql|genome_statistics_value_longer');
