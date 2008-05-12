# patch_49_50_b.sql
#
# title: change defaults in coord_system
#
# description:
# Make version column in coord_system default to null

ALTER TABLE coord_system CHANGE COLUMN version version VARCHAR(40) DEFAULT NULL;

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_49_50_b.sql|coord_system_version_default');


