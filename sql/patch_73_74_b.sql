# patch_73_74_b.sql
#
# title: dnac removal
#
# description:
# Removal dnac table which is not used any more

DROP TABLE dnac;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_73_74_b.sql|remove_dnac');

 
