# patch_63_64_c.sql
#
# Title:
#   is_ref to alt_allele 
#
# Description:
#    Add is_ref to alt_allele to specify which is the reference al


# Add the new column.
ALTER TABLE alt_allele
  ADD COLUMN is_ref BOOLEAN NOT NULL DEFAULT 0
    AFTER gene_id;

# Insert patch identifier.
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch',
    'patch_63_64_c.sql|is_ref_added_to_alt_allele');
