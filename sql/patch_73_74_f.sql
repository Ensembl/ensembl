# patch_73_74_f.sql
#
# title: Pair_dna_align_feature removal
#
# description:
# Removal of the pair_dna_align_feature_id column in the dna_align_feature table, as it is not used any more

ALTER TABLE dna_align_feature DROP INDEX pair_idx;
ALTER TABLE dna_align_feature DROP COLUMN pair_dna_align_feature_id ;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_73_74_f.sql|remove_pair_dna_align');

 
