# patch_50_51_c.sql
#
# title: Change indices on meta_coord
#
# Description:
#   Swap order of table_name,coord_system_id ID in meta_coord unique
#   index.

DROP INDEX table_name ON meta_coord;

CREATE UNIQUE INDEX cs_table_name_idx
  ON meta_coord (coord_system_id, table_name);

ANALYZE TABLE meta_coord;

# Patch identifier
INSERT INTO meta (meta_key, meta_value)
  VALUES ('patch', 'patch_50_51_c.sql|meta_coord_indices');
