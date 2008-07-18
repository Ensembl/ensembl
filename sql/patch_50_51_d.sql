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