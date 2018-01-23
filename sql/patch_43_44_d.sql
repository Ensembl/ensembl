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

# patch_43_44_d
#
# title: Translation stable_id unique
#
# description:
# Remove UNIQUE constraint from translation stable ID to make it the same as other stable_id tables.

ALTER TABLE translation_stable_id DROP INDEX stable_id;
ALTER TABLE translation_stable_id ADD INDEX stable_id_idx (stable_id, version);

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_43_44_d.sql|translation_stable_id_unique');

