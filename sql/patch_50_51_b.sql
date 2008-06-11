# patch_50_51_b.sql
#
# Title: make database multi-species capable
#
# Description:
#   Add the species_meta table and move all species-specific meta data
#   from the meta table into this new table.  Also add a species_id
#   column to the coord_system table.

-- Create the new table species_meta with the meta table as the
-- reference
CREATE TABLE species_meta LIKE meta;

-- Drop the indexes inherited from meta
ALTER TABLE species_meta DROP INDEX key_value;
ALTER TABLE species_meta DROP INDEX meta_key_index;
ALTER TABLE species_meta DROP INDEX meta_value_index;

-- Change the name of meta_id into species_meta_id on species_meta
ALTER TABLE species_meta CHANGE COLUMN meta_id
  species_meta_id INT(10) UNSIGNED NOT NULL AUTO_INCREMENT;

-- Add the new species_id column after species_meta_id
ALTER TABLE species_meta ADD COLUMN
  species_id INT(10) UNSIGNED NOT NULL DEFAULT 1    -- Default species_id is 1
  AFTER species_meta_id;

-- Add indexes to species_meta
ALTER TABLE species_meta ADD INDEX meta_key_idx (meta_key, species_id);

-- Add the species_id column to coord_system table
ALTER TABLE coord_system ADD COLUMN
  species_id INT(10) UNSIGNED NOT NULL DEFAULT 1    -- Default species_id is 1
  AFTER coord_system_id;

-- Drop the indexes from coord_system
ALTER TABLE coord_system DROP INDEX rank;
ALTER TABLE coord_system DROP INDEX name;

-- Add new indexes to coord_system
ALTER TABLE coord_system ADD UNIQUE INDEX rank_idx (rank, species_id);
ALTER TABLE coord_system ADD UNIQUE INDEX name_idx (name, version, species_id);
ALTER TABLE coord_system
  ADD UNIQUE INDEX species_cs_idx (species_id, coord_system_id);

-- Move data from the meta table into the species_meta table
INSERT INTO species_meta (species_id, meta_key, meta_value)
  SELECT    1, m.meta_key, m.meta_value             -- Default species_id is 1
  FROM      meta m
  WHERE     m.meta_key != 'patch'
    AND     m.meta_key != 'schema_version'
  ORDER BY  m.meta_key, m.meta_id;

DELETE FROM meta
  WHERE     meta_key != 'patch'
    AND     meta_key != 'schema_version';

-- Optimize the modified tables
OPTIMIZE TABLE meta;
OPTIMIZE TABLE species_meta;
OPTIMIZE TABLE coord_system;

# patch identifier
INSERT INTO meta (meta_key, meta_value)
VALUES ('patch', 'patch_50_51_b.sql|multispecies');
