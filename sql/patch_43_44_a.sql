# patch_43_44_a
#
# title: Key column types
#
# description:
# Change the types of some columns to make them more appropriate to the values they store.

ALTER TABLE analysis CHANGE COLUMN analysis_id analysis_id SMALLINT UNSIGNED NOT NULL AUTO_INCREMENT;
ALTER TABLE oligo_feature CHANGE COLUMN analysis_id analysis_id SMALLINT UNSIGNED NOT NULL;
ALTER TABLE analysis_description CHANGE COLUMN analysis_id analysis_id SMALLINT UNSIGNED NOT NULL;
ALTER TABLE simple_feature CHANGE COLUMN analysis_id analysis_id SMALLINT UNSIGNED NOT NULL;
ALTER TABLE protein_align_feature CHANGE COLUMN analysis_id analysis_id SMALLINT UNSIGNED NOT NULL;
ALTER TABLE dna_align_feature CHANGE COLUMN analysis_id analysis_id SMALLINT UNSIGNED NOT NULL;
ALTER TABLE repeat_feature CHANGE COLUMN analysis_id analysis_id SMALLINT UNSIGNED NOT NULL;
ALTER TABLE gene CHANGE COLUMN analysis_id analysis_id SMALLINT UNSIGNED NOT NULL;
ALTER TABLE transcript CHANGE COLUMN analysis_id analysis_id SMALLINT UNSIGNED NOT NULL;
ALTER TABLE protein_feature CHANGE COLUMN analysis_id analysis_id SMALLINT UNSIGNED NOT NULL;
ALTER TABLE identity_xref CHANGE COLUMN analysis_id analysis_id SMALLINT UNSIGNED NOT NULL;
ALTER TABLE prediction_transcript CHANGE COLUMN analysis_id analysis_id SMALLINT UNSIGNED NOT NULL;
ALTER TABLE marker_feature CHANGE COLUMN analysis_id analysis_id SMALLINT UNSIGNED NOT NULL;
ALTER TABLE qtl_feature CHANGE COLUMN analysis_id analysis_id SMALLINT UNSIGNED NOT NULL;
ALTER TABLE density_type CHANGE COLUMN analysis_id analysis_id SMALLINT UNSIGNED NOT NULL;
ALTER TABLE regulatory_feature CHANGE COLUMN analysis_id analysis_id SMALLINT UNSIGNED NOT NULL;
ALTER TABLE regulatory_search_region CHANGE COLUMN analysis_id analysis_id SMALLINT UNSIGNED NOT NULL;
ALTER TABLE unmapped_object CHANGE COLUMN analysis_id analysis_id SMALLINT UNSIGNED NOT NULL;
ALTER TABLE ditag_feature CHANGE COLUMN analysis_id analysis_id SMALLINT UNSIGNED NOT NULL;

ALTER TABLE external_db CHANGE COLUMN external_db_id external_db_id SMALLINT UNSIGNED NOT NULL;
ALTER TABLE xref CHANGE COLUMN external_db_id external_db_id SMALLINT UNSIGNED NOT NULL;
ALTER TABLE unmapped_object CHANGE COLUMN external_db_id external_db_id SMALLINT UNSIGNED NOT NULL;

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_43_44_a.sql|rationalise_key_columns');

