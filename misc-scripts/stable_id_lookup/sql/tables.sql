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

CREATE TABLE species (
  species_id		INTEGER UNSIGNED NOT NULL AUTO_INCREMENT,
  name         	     	VARCHAR(255) NOT NULL,
  taxonomy_id        	INTEGER UNSIGNED,

  PRIMARY KEY (species_id),
  UNIQUE INDEX name_idx (name)
) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


CREATE TABLE stable_id_lookup (
 stable_id   	  VARCHAR(128) NOT NULL,	      
 species_id	  INTEGER UNSIGNED NOT NULL,
 db_type          VARCHAR(255) NOT NULL,
 object_type   	  VARCHAR(255) NOT NULL
) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


CREATE TABLE archive_id_lookup (
 archive_id        VARCHAR(128) NOT NULL,
 species_id       INTEGER UNSIGNED NOT NULL,
 db_type          VARCHAR(255) NOT NULL,
 object_type      VARCHAR(255) NOT NULL,

 UNIQUE INDEX archive_id_lookup_idx (archive_id,species_id,db_type,object_type),
 KEY archive_id_db_type (archive_id,db_type,object_type),
 KEY archive_id_object_type (archive_id,object_type)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


CREATE TABLE meta (

  meta_id                     INT NOT NULL AUTO_INCREMENT,
  species_id                  INT UNSIGNED DEFAULT 1,
  meta_key                    VARCHAR(40) NOT NULL,
  meta_value                  VARCHAR(255) BINARY,

  PRIMARY   KEY (meta_id),
  UNIQUE    KEY species_key_value_idx (species_id, meta_key, meta_value),
            KEY species_value_idx (species_id, meta_value)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

