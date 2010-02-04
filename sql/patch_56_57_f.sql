# patch_56_57_f.sql
#
# Title: Expand width of simple_feature.display_label
#
# Description:
# The display_label field of the simple_feature table needs to be
# expanded from 40 characters to be able to hold some of the data
# imported from Encode for the human database.

ALTER TABLE simple_feature
MODIFY COLUMN display_label VARCHAR(255) NOT NULL;

OPTIMIZE TABLE simple_feature;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
VALUES (NULL, 'patch', 'patch_56_57_f.sql|simple_feature.display_label');
