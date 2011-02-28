# patch_61_62_c.sql
#
# Title: Index for db_name and db_release.
#
# Description:
# Add unique index to db_name and db_release fields in external_db table.

CREATE UNIQUE INDEX db_name_db_release_idx ON external_db(db_name,db_release);

# Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_61_62_c.sql|db_name_db_release_idx');
