# patch_73_74_f.sql
#
# title: Canonical_annotation removal
#
# description:
# Removal of the canonical_annotation column in the gene table, as it is not used any more

ALTER TABLE gene DROP COLUMN canonical_annotation ;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_73_74_e.sql|remove_canonical_annotation');

 
