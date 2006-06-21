# patch_39_40_b
#
# title: xref unique constraint
#
# description:
# Add info_type and info_text columns to the id_index UNIQUE KEY in the xref table.

ALTER TABLE xref DROP INDEX id_index;
ALTER TABLE xref ADD UNIQUE id_index (dbprimary_acc, external_db_id, info_type, info_text);

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_39_40_b.sql|xref_unique_constraint');

