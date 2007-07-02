# patch_45_46_d.sql
#
# title: Meta unique key
#
# description:
# Add unique key to meta_key/meta_value columns in meta table.

ALTER TABLE meta ADD UNIQUE key_value (meta_key, meta_value);

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_45_46_d.sql|meta_unique_key');



