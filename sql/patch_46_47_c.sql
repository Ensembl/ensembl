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

# patch_46_47_c.sql
#
# title: extend external_db db_release and db_name
#
# description:
# Make db_release column of external_db VARCHAR(255)

ALTER TABLE external_db CHANGE COLUMN db_release db_release VARCHAR(255);
ALTER TABLE external_db CHANGE COLUMN db_name db_name VARCHAR(28) NOT NULL;

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_46_47_c.sql|extend_db_release');


