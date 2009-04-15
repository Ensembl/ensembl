# patch_50_51_e.sql
#
# title: Add external_data column to protein_feature and dna_align_feature.
#
# description:
# Add external_data column to protein_feature and dna_align_feature tables, primarily for storage of user-uploaded data.

ALTER TABLE protein_feature ADD COLUMN external_data TEXT;

ALTER TABLE dna_align_feature ADD COLUMN external_data TEXT;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_50_51_e.sql|feature_external_data');


