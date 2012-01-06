# patch_65_66_e.sql
#
# Title: Make xref.external_db_id INTEGER UNSIGNED NOT NULL
#
# Description:
# Add NOT NULL to definition of external_db_id
ALTER TABLE xref                    MODIFY external_db_id INTEGER UNSIGNED NOT NULL;

# Patch identifier:
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_65_66_e.sql|fix_external_db_id_in_xref');
