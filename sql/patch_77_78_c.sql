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

#
# Title: patch_77_78_c.sql - Unmapped_reason_id
#
# Description: Change unmapped_reason_id from smallint to int
#   

ALTER TABLE unmapped_reason MODIFY COLUMN unmapped_reason_id int(10) unsigned NOT NULL AUTO_INCREMENT;
ALTER TABLE unmapped_object MODIFY COLUMN unmapped_reason_id int(10) unsigned NOT NULL;

-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_77_78_c.sql|Change unmapped_reason_id from smallint to int');
