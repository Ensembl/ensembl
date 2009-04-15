# patch_50_51_b.sql
#
# title: protein_feature hit_name
#
# description:
# Rename hit_id column in protein_feature to hit_name

ALTER TABLE protein_feature DROP INDEX hid_index;

ALTER TABLE protein_feature CHANGE COLUMN hit_id hit_name VARCHAR(40) NOT NULL;

ALTER TABLE protein_feature ADD INDEX hitname_idx (hit_name);

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_50_51_b.sql|protein_feature_hit_name');


