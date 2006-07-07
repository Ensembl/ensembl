# patch_39_40_e
#
# title: marker_synonym name
#
# description:
# this patch widens the name column in marker_synonym to 50 characters.

ALTER TABLE marker_synonym CHANGE name name VARCHAR(50);

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_39_40_j.sql|marker_synonym_name');

