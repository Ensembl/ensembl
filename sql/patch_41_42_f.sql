# patch_40_41_f
#
# title: Analysis_description.web_data
#
# description: Add web_data column to analysis_description.

ALTER TABLE analysis_description ADD COLUMN web_data TEXT;

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_41_42_f.sql|analysis_description_web_data');
