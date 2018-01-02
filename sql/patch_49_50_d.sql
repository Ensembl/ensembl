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

# patch_49_50_d.sql
#
# title: Change indices on seq_region
#
# Description: Swap order of name,coord_system ID in seq_region unique index, and add new index on coord_system ID. Remove redundant name index.

DROP INDEX name_idx         ON seq_region;
DROP INDEX coord_system_id  ON seq_region;

CREATE UNIQUE INDEX name_cs_idx ON seq_region (name, coord_system_id);

CREATE INDEX cs_idx ON seq_region (coord_system_id);

ANALYZE TABLE seq_region;

# Patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_49_50_d.sql|seq_region_indices');



