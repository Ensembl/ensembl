# patch_46_47_b.sql
#
# title: new align columns
#
# description:
# Add external_db_id and hcoverage columns to dna_align_feature and protein_align_feature

ALTER TABLE dna_align_feature ADD COLUMN external_db_id SMALLINT UNSIGNED;
ALTER TABLE dna_align_feature ADD COLUMN hcoverage DOUBLE;

ALTER TABLE dna_align_feature ADD KEY external_db_idx (external_db_id);

ALTER TABLE protein_align_feature ADD COLUMN external_db_id SMALLINT UNSIGNED;
ALTER TABLE protein_align_feature ADD COLUMN hcoverage DOUBLE;

ALTER TABLE protein_align_feature ADD KEY external_db_idx (external_db_id);

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_46_47_b.sql|new_align_columns');


