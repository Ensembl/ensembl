# patch_50_51_g.sql
#
# title: Protein feature score.
#
# description:
# Allow the score column in protein_feature to be null.

ALTER TABLE `protein_feature` MODIFY COLUMN `score` DOUBLE;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_50_51_g.sql|protein_feature_score');


