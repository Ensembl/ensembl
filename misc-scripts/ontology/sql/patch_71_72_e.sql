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

-- patch_71_72_e.sql
--
-- Title: Add is_obsolete
--
-- Description:
--   Adds the is_obsolete flag to the term table

ALTER TABLE TERM
ADD COLUMN is_obsolete INT NOT NULL DEFAULT 0;

-- Patch identifier
INSERT INTO meta (meta_key, meta_value)
  VALUES ('patch', 'patch_71_72_e.sql|is_obsolete');


