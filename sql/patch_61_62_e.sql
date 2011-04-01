# patch_61_62_e.sql
#
# Title: Index for seq_region_synonym.seq_region_id
#
# Description:
# Add an index to the seq_region_id column in seq_region_synonym

CREATE INDEX seq_region_idx ON seq_region_synonym(seq_region_id);

# Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_61_62_e.sql|seq_region_synonym_seq_region_idx');
