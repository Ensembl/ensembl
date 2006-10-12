# patch_41_42_c
#
# title: unique analysis id in analysis_description
#
# description:
# Add a UNIQUE constraint to the analysis_description.analysis_id

ALTER TABLE analysis_description DROP INDEX analysis_idx;
ALTER TABLE analysis_description ADD UNIQUE analysis_idx (analysis_id);

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_41_42_c.sql|analysis_description_unique');