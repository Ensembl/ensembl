-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2018] EMBL-European Bioinformatics Institute
-- 
-- Licensed under the Apache License, Version 2.0 (the "License");
-- you may not use this file except in compliance with the License.
-- You may obtain a copy of the License at
-- 
--      http://www.apache.org/licenses/LICENSE-2.0
-- 
-- Unless required by applicable law or agreed to in writing, software
-- distributed under the License is distributed on an "AS IS" BASIS,
-- WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
-- See the License for the specific language governing permissions and
-- limitations under the License.

# patch_71_72_b.sql
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
  KEY associated_group_idx (associated_group_id),
  UNIQUE KEY object_associated_source_type_idx (object_xref_id, xref_id, source_xref_id, condition_type, associated_group_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


CREATE TABLE associated_group (

  associated_group_id            INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  description                    VARCHAR(128) DEFAULT NULL,

  PRIMARY KEY (associated_group_id)
) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

# Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_71_72_b.sql|associated_xref');

