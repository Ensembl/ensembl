# patch_50_51_c.sql
#
# title: meta_coord index
#
# description:
# Swap order of columns in meta_coord table index for better performance.

DROP INDEX table_name ON meta_coord;

CREATE UNIQUE INDEX cs_table_name_idx ON meta_coord (coord_system_id, table_name);

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_50_51_c.sql|meta_coord_index');


