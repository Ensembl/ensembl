# patch_51_52_c.sql
#
# title: Add dna_align_feature.pair_dna_align_feature_id
#
# description:
# Add a pair_dna_align_feature_id to allow support for paired reads

ALTER TABLE dna_align_feature ADD COLUMN pair_dna_align_feature_id INT(10) UNSIGNED;

ALTER TABLE dna_align_feature ADD INDEX pair_idx (pair_dna_align_feature_id);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_51_52_c.sql|pair_dna_align_feature_id');


