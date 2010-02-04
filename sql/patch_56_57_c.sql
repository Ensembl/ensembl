# patch_56_57_c.sql
#
# title: Modify the external_db.type enumeration
#
# description:
# Modify the external_db.type enumeration so that it also includes 'ALT_GENE'.

ALTER TABLE external_db MODIFY type ENUM('ARRAY', 'ALT_TRANS', 'ALT_GENE', 'MISC', 'LIT', 'PRIMARY_DB_SYNONYM', 'ENSEMBL');

# patch identifier

INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_56_57_c.sql|external_db_type_enum');
