# patch_43_44_d
#
# title: Translation stable_id unique
#
# description:
# Remove UNIQUE constraint from translation stable ID to make it the same as other stable_id tables.

ALTER TABLE translation_stable_id DROP INDEX stable_id;
ALTER TABLE translation_stable_id ADD INDEX stable_id_idx (stable_id, version);

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_43_44_d.sql|translation_stable_id_unique');

