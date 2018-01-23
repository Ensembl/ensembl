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

-- patch_71_72_d.sql
--
-- Title: Fix patch versions
--
-- Description:
--   Fixes the existing patch meta items as their versioning was not correct

update meta set meta_value = 'patch_71_72_b.sql|alt_id table'
  where meta_key = 'patch'
  and meta_value = 'patch_71_72b.sql|alt_id table';

update meta set meta_value = 'patch_71_72_c.sql|schema_version'
  where meta_key = 'patch'
  and meta_value = 'patch_71_72c.sql|schema_version';

-- Patch identifier
INSERT INTO meta (meta_key, meta_value)
  VALUES ('patch', 'patch_71_72_d.sql|patch_version_fix');


