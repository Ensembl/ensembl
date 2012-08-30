# patch_67_68_b.sql
#
# Title: Create associated_xref table
#
# Description:
#  Create table associated_xref for associating object xrefs with an associated
#  annotation (eg Gene Ontology Annotation Extensions) given a source xref and 
#  condition.

CREATE TABLE associated_xref (
  
  associated_xref_id             INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,  
  object_xref_id                 INT(10) UNSIGNED DEFAULT '0' NOT NULL,
  xref_id                        INT(10) UNSIGNED DEFAULT '0' NOT NULL,
  source_xref_id                 INT(10) UNSIGNED DEFAULT NULL,
  condition_type                 VARCHAR(128) DEFAULT NULL,
  associated_group_id            INT(10) UNSIGNED DEFAULT NULL,
  rank                           INT(10) UNSIGNED DEFAULT '0',
  
  PRIMARY KEY (associated_xref_id),
  KEY associated_source_idx (source_xref_id),
  KEY associated_object_idx (object_xref_id),
  KEY associated_idx (xref_id),
  FOREIGN KEY associated_group_idx (associated_group_id) REFERENCES associated_group (associated_group_id),
  UNIQUE KEY object_associated_source_type_idx (object_xref_id, xref_id, source_xref_id, condition_type, associated_group_id)
  
) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


CREATE TABLE associated_group (
  
  associated_group_id            INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,  
  description                    VARCHAR(128) DEFAULT NULL,
  
  PRIMARY KEY (associated_group_id)
) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

# Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_67_68_x.sql|associated_xref');
