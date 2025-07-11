-- See the NOTICE file distributed with this work for additional information
-- regarding copyright ownership.
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

# patch_45_46_c.sql
#
# title: unmapped_object.external_db_id
#
# description:
# Make external_db column in unmapped_object be SMALLINT UNSIGNED explicitly.

ALTER TABLE unmapped_object CHANGE COLUMN external_db_id external_db_id SMALLINT UNSIGNED;

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_45_46_c.sql|unmapped_object.external_db_id');



