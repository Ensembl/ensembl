# patch_65_66_d.sql

# Title: Index for ontology_xref.object_xref_id
#
# Description:
# Add an index to the object_xref_id column in ontology_xref

CREATE INDEX object_idx ON ontology_xref(object_xref_id);

# Patch identifier:
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_65_66_d.sql|add_index_to_ontology_xref_table');
