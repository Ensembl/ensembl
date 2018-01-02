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

CREATE TABLE external_db_type (
  external_db_id	INTEGER UNSIGNED NOT NULL AUTO_INCREMENT,
  db_name       	VARCHAR(100) NOT NULL,
  db_display_name	VARCHAR(255),
  db_type_id		INTEGER UNSIGNED,

  PRIMARY KEY (external_db_id),
  KEY type_idx (db_type_id, db_name)
) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


CREATE TABLE object_xref_linkage (
 external_db_id   	  INTEGER UNSIGNED NOT NULL,
 ensembl_object_type	  VARCHAR(100),
 link_type_id	  	  INTEGER UNSIGNED NOT NULL,
 linked_external_db_id    INTEGER UNSIGNED,
 linked_node_text	  VARCHAR(255),

 KEY xref_origin_idx (external_db_id, ensembl_object_type)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

CREATE TABLE object_external_db_node (
 ensembl_object_type      VARCHAR(100) NOT NULL,
 external_db_id		  INTEGER UNSIGNED NOT NULL,
 mindmap_tag_id		  VARCHAR(100) NOT NULL,

 UNIQUE KEY (ensembl_object_type, external_db_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


CREATE TABLE db_type (
  db_type_id    INTEGER UNSIGNED NOT NULL AUTO_INCREMENT,
  db_type	VARCHAR(255) NOT NULL,

  PRIMARY KEY (db_type_id)
) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

CREATE TABLE link_type (
  link_type_id	       INTEGER UNSIGNED NOT NULL AUTO_INCREMENT,
  link_type	       VARCHAR(100) NOT NULL,
  link_description     VARCHAR(255) NOT NULL,

  PRIMARY KEY (link_type_id),
  KEY link_type_idx (link_type)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

CREATE TABLE object_distance (
  from_object	       VARCHAR(100) NOT NULL,
  to_object	       VARCHAR(100) NOT NULL,
  distance	       TINYINT(1) UNSIGNED NOT NULL,

  KEY (from_object,distance)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

LOCK TABLES `object_distance` WRITE;
INSERT INTO `object_distance` VALUES ('Gene','Gene',0), ('Gene','Transcript',1), ('Gene', 'Translation', 2), ('Transcript','Gene',1), ('Transcript','Transcript',0), ('Transcript', 'Translation', 1), ('Translation', 'Gene', 2), ('Translation', 'Transcript',1), ('Translation', 'Translation',0);
UNLOCK TABLES;

LOCK TABLES `db_type` WRITE;
INSERT INTO `db_type` VALUES (1,'disease related'),(2,'expression'),(3,'integrated information'),(4,'naming'),(5,'products'),(6,'sequence/annotation'),(7,'classification'),(8,'function/location'),(9,'structure'),(10,'other resources');
UNLOCK TABLES;

LOCK TABLES `link_type` WRITE;
INSERT INTO `link_type` VALUES (1,'DIRECT', 'DIRECT'),(2,'INFERRED PAIR', 'INFERRED_PAIR'),(3,'DEPENDENT', 'DEPENDENT ON'),(4,'SEQUENCE_MATCH', 'SEQUENCE MATCH'),(5,'COORDINATE_OVERLAP', 'COORDINATE OVERLAP'),(6,'GENERATED_FROM', 'GENERATED FROM'),(7,'PROJECTION','PROJECTION'),(8,'PROTEIN_FEATURES', 'VIA PROTEIN FEATURES');
UNLOCK TABLES;



