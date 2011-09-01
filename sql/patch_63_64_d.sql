# patch_63_64_d.sql
#
# Title:
#   change to linkage_type column in ontology_xref 
#
# Description:
#    Change the column definition to VARCHAR(3) instead of providing
#    a list of allowable values. Values will be tested in healthchecks.

ALTER TABLE ontology_xref DROP INDEX object_source_type_idx;

ALTER TABLE ontology_xref
  MODIFY linkage_type VARCHAR(3) DEFAULT NULL;

ALTER TABLE ontology_xref 
  ADD UNIQUE INDEX object_source_type_idx (object_xref_id, source_xref_id, linkage_type);

# Insert patch identifier.
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch',
    'patch_63_64_d.sql|linkage_type change in ontology_xref');
