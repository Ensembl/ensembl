# patch_74_75_e.sql
#
# title: Add unique constraint on attrib tables
#
# description:
# For all attrib related tables, attrib_type_id and value should be unique for a given object

ALTER TABLE seq_region_attrib ADD UNIQUE KEY region_attribx (seq_region_id, attrib_type_id, value(40));
ALTER TABLE gene_attrib ADD UNIQUE KEY gene_attribx (gene_id, attrib_type_id, value(40));
ALTER TABLE transcript_attrib ADD UNIQUE KEY transcript_attribx (transcript_id, attrib_type_id, value(40));
ALTER TABLE translation_attrib ADD UNIQUE KEY translation_attribx (translation_id, attrib_type_id, value(40));
ALTER TABLE misc_attrib ADD UNIQUE KEY misc_attribx (misc_feature_id, attrib_type_id, value(40));


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_74_75_e.sql|unique_attrib_key');

 
