# patch_43_44_f
#
# title: Add PRIMARY_DB_SYN to type
#
# description:
# Add an extra value to the type ENUM for external_db

ALTER TABLE external_db CHANGE COLUMN type type ENUM('ARRAY', 'ALT_TRANS', 'MISC', 'LIT', 'PRIMARY_DB_SYNONYM');

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_43_44_f.sql|external_db_type_syn');

