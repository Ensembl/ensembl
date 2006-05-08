# patch_38_39_d
#
# title: schema version
#
# description:
# this patch updates the schema version

# update schema version
UPDATE meta set meta_value = 39 where meta_key = 'schema_version';

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_38_39_d.sql|schema_version');

