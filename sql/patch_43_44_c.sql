# patch_43_44_c
#
# title: Add type column to external_db
#
# description:
# Add a column for the type of the xref (array or whatever) to external_db

ALTER TABLE external_db ADD COLUMN type ENUM('ARRAY', 'ALT_TRANS', 'MISC', 'LIT');

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_43_44_c.sql|external_db_type');

