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

# patch_44_45_c.sql
#
# title: db_release not null
#
# description: 
# Remove NOT NULL constraint on external_db.db_release

ALTER TABLE external_db CHANGE COLUMN db_release db_release VARCHAR(40);

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_44_45_c.sql|db_release_not_null');


