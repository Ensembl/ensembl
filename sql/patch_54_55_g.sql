# patch_54_55_g.sql
#
# title: Change analysis_description.display_label to "NOT NULL".
#
# description:
# Display labels are generally supposed to be non-NULL and non-empty.
# The display_label field in the analysis_description table allows for
# NULLs.  This patch fixes this.

ALTER TABLE analysis_description
  MODIFY COLUMN display_label VARCHAR(255) NOT NULL;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
VALUES (NULL, 'patch',
  'patch_54_55_g.sql|analysis_description.display_label_NOT_NULL');
