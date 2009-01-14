# patch_52_53_b.sql
#
# title: external_db type enum
#
# description:
# Add the enum value 'ENSEMBL' to the external_db.type column (primarily for eFG)

ALTER TABLE external_db MODIFY type ENUM('ARRAY', 'ALT_TRANS', 'MISC', 'LIT', 'PRIMARY_DB_SYNONYM', 'ENSEMBL');

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_52_53_b.sql|external_db_type_enum');


