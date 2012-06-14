# patch_67_68_b.sql
#
# Title: Create associated_xref table
#
# Description:
#  Create table associated_xref for associating object xrefs with an associated
#  annotation (eg Gene Ontology Annotation Extensions) given a source xref and 
#  condition.

CREATE TABLE associated_xref (
  
  object_xref_id                 INT(10) UNSIGNED DEFAULT '0' NOT NULL,
  associated_xref_id             INT(10) UNSIGNED DEFAULT NULL,
  source_xref_id                 INT(10) UNSIGNED DEFAULT NULL,
  condition_type                 VARCHAR(128) DEFAULT NULL,
  linked_associated_xref_id      INT(10) UNSIGNED DEFAULT '0' NOT NULL,
  
  KEY associated_source_idx (source_xref_id),
  KEY associated_object_idx (object_xref_id),
  KEY associated_idx (associated_xref_id),
  KEY linked_associated_idx (linked_associated_xref_id),
  UNIQUE KEY object_associated_source_type_idx (object_xref_id, associated_xref_id, source_xref_id, condition_type, linked_associated_xref_id)
  
) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

# Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_67_68_x.sql|associated_xref');
