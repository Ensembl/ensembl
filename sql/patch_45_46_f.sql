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

# patch_45_46_f.sql
#
# title stable_id_event.uni_idx
#
# description:
#   Drop old_version and new_version from uni_idx on stable_id_event

ALTER TABLE stable_id_event DROP KEY uni_idx;

ALTER TABLE stable_id_event ADD UNIQUE KEY uni_idx
  ( mapping_session_id, old_stable_id, new_stable_id, type );

# patch identifier
INSERT INTO meta (meta_key, meta_value) 
  VALUES ('patch', 'patch_45_46_f.sql|stable_id_event.uni_idx');



