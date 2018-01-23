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

# patch_65_66_b.sql
#
# Title: Make external_db.external_db_id AUTO_INCREMENT and INTEGER UNSIGNED.
#
# Description:
# We're using too high values in external_db.external_db_id for the
# current SMALLINT, and with the web interface we're using internally
# to add new entries, we also need this field to be AUTO_INCREMENT.

ALTER TABLE external_db
  MODIFY external_db_id INTEGER UNSIGNED NOT NULL AUTO_INCREMENT;

# Also modify this field in the other tables that uses it as a foreign key:
ALTER TABLE dna_align_feature       MODIFY external_db_id INTEGER UNSIGNED;
ALTER TABLE protein_align_feature   MODIFY external_db_id INTEGER UNSIGNED;
ALTER TABLE seq_region_synonym      MODIFY external_db_id INTEGER UNSIGNED;
ALTER TABLE unmapped_object         MODIFY external_db_id INTEGER UNSIGNED;
ALTER TABLE xref                    MODIFY external_db_id INTEGER UNSIGNED;

# Patch identifier:
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_65_66_b.sql|fix_external_db_id');
