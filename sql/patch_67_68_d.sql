# patch_67_68_d.sql
#
# Title: Altering intron_supporting_evidence table
#
# Description: Adding is_splice_canonical which lets us define if a splice 
#              junction can be considered canonical and adding missing indexes
# 

ALTER TABLE intron_supporting_evidence 
ADD COLUMN is_splice_canonical BOOLEAN NOT NULL DEFAULT 0;

ALTER TABLE intron_supporting_evidence ADD KEY seq_region_idx (seq_region_id, seq_region_start);

# Patch identifier:
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_67_68_d.sql|add_is_splice_canonical_and_seq_index');