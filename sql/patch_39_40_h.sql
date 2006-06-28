# patch_39_40_h
#
# title: oligo_feature analysis id column type change
#
# description:
# Change oligo_feature.analysis_id to be int(10) unsigned. Should have been part of patch_39_40_a.sql

ALTER TABLE oligo_feature CHANGE COLUMN analysis_id analysis_id INT(10) UNSIGNED NOT NULL;

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_39_40_h.sql|oligo_feature_analysis_id_type');

