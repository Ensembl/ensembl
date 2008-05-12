# patch_49_50_d.sql
#
# title: Change indices on seq_region
#
# Description: Swap order of name,coord_system ID in seq_region unique index, and add new index on coord_system ID. Remove redundant name index.

DROP INDEX name_idx         ON seq_region;
DROP INDEX coord_system_id  ON seq_region;

CREATE UNIQUE INDEX name_cs_idx ON seq_region (name, coord_system_id);

CREATE INDEX cs_idx ON seq_region (coord_system_id);

ANALYZE TABLE seq_region;

# Patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_49_50_d.sql|seq_region_indices');



