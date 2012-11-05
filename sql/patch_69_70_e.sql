# patch_69_70_e.sql
#
# Title: Add hit_description to protein_feature
#
# Description: Column allowing an optional description to be added to a protein_feature

ALTER TABLE protein_feature ADD COLUMN hit_description TEXT;

# Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_69_70_e.sql|protein_feature_hit_description');


