# patch_73_74_d.sql
#
# title: Qtl removal
#
# description:
# Removal of the qtl tables (qtl, qtl_feature, qtl_synonym) which are not used any more

DROP TABLE qtl;
DROP TABLE qtl_feature;
DROP TABLE qtl_synonym;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_73_74_d.sql|remove_qtl');

 
