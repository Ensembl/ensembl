# patch_61_62_b.sql
#
# Title: 'Synonym' field extension.
#
# Description:
#   Extend field 'synonym' in 'external_synonym' table to 100 chars
#   to support some very long synonyms in D. melanogaster.

ALTER table external_synonym MODIFY synonym VARCHAR(100) NOT NULL;

# Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_61_62_b.sql|synonym_field_extension');
