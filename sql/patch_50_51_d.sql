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

# patch_50_51_d.sql
#
# Title: make database multi-species capable
#
# Description:
#   Add a species_id column to the meta and coord_system table and make
#   new indexes on these tables.

-- Add the new species_id column after meta_id
ALTER TABLE meta ADD COLUMN
 species_id INT UNSIGNED DEFAULT 1 -- Default species_id is 1
                                   -- NULL means "not species specific"
 AFTER meta_id;

-- Redo the indexes on the meta table
ALTER TABLE meta DROP INDEX key_value;
ALTER TABLE meta DROP INDEX meta_key_index;
ALTER TABLE meta DROP INDEX meta_value_index;

ALTER TABLE meta
 ADD UNIQUE INDEX species_key_value_idx (species_id, meta_key, meta_value);
ALTER TABLE meta
 ADD INDEX species_value_idx (species_id, meta_value);

-- Add the species_id column to coord_system table
-- Default species_id is 1
ALTER TABLE coord_system ADD COLUMN species_id INT(10) UNSIGNED NOT NULL DEFAULT 1 AFTER coord_system_id;

-- Drop the indexes from coord_system
ALTER TABLE coord_system DROP INDEX rank;
ALTER TABLE coord_system DROP INDEX name;

-- Add new indexes to coord_system
ALTER TABLE coord_system ADD UNIQUE INDEX rank_idx (rank, species_id);
ALTER TABLE coord_system ADD UNIQUE INDEX name_idx (name, version, species_id);
ALTER TABLE coord_system ADD        INDEX species_idx (species_id);

-- Optimize the modified tables
OPTIMIZE TABLE meta;
OPTIMIZE TABLE coord_system;

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_50_51_d.sql|multispecies');
