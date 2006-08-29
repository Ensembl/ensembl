# patch_40_41_a
#
# title: analysis_description.displayable
#
# description:
# Add displayable column to analysis_description

ALTER TABLE analysis_description ADD COLUMN displayable BOOLEAN DEFAULT 1 NOT NULL;

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_40_41_a.sql|analysis_description_displayable');

