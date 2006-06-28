# patch_39_40_i
#
# title: schema version
#
# description:
# this patch updates the schema version

# update schema version
UPDATE meta set meta_value = 40 where meta_key = 'schema_version';

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_39_40_i.sql|schema_version');

