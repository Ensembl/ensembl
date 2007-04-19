# patch_44_45_c.sql
#
# title: db_release not null
#
# description: 
# Remove NOT NULL constraint on external_db.db_release

ALTER TABLE external_db CHANGE COLUMN db_release db_release VARCHAR(40);

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_44_45_c.sql|db_release_not_null');


