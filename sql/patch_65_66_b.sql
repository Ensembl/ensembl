# patch_65_66_b.sql
#
# Title: Make external_db.external_db_id AUTO_INCREMENT and INTEGER UNSIGNED.
#
# Description:
# We're using too high values in external_db.external_db_id for the
# current SMALLINT, and with the web interface we're using internally
# to add new entries, we also need this field to be AUTO_INCREMENT.

ALTER TABLE external_db
  MODIFY external_db_id INTEGER UNSIGNED NOT NULL AUTO_INCREMENT;

# Also modify this field in the other tables that uses it as a foreign key:
ALTER TABLE dna_align_feature       MODIFY external_db_id INTEGER UNSIGNED;
ALTER TABLE protein_align_feature   MODIFY external_db_id INTEGER UNSIGNED;
ALTER TABLE seq_region_synonym      MODIFY external_db_id INTEGER UNSIGNED;
ALTER TABLE unmapped_object         MODIFY external_db_id INTEGER UNSIGNED;
ALTER TABLE xref                    MODIFY external_db_id INTEGER UNSIGNED;

# Patch identifier:
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_65_66_b.sql|fix_external_db_id');
