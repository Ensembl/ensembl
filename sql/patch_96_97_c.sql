-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2020] EMBL-European Bioinformatics Institute
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

# patch_96_97_c.sql
#
# Title: Create tables for rnaproduct data
#
# Description:
#   These tables will be used to store information about mature RNA products, e.g. microRNA

CREATE TABLE rnaproduct (

  rnaproduct_id               INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  rnaproduct_type_id          SMALLINT(5) UNSIGNED NOT NULL,
  transcript_id               INT(10) UNSIGNED NOT NULL,
  seq_start                   INT(10) NOT NULL,       # relative to transcript start
  start_exon_id               INT(10) UNSIGNED,
  seq_end                     INT(10) NOT NULL,       # relative to transcript start
  end_exon_id                 INT(10) UNSIGNED,
  stable_id                   VARCHAR(128) DEFAULT NULL,
  version                     SMALLINT UNSIGNED DEFAULT NULL,
  created_date                DATETIME DEFAULT NULL,
  modified_date               DATETIME DEFAULT NULL,

  PRIMARY KEY (rnaproduct_id),
  KEY transcript_idx (transcript_id),
  KEY stable_id_idx (stable_id, version)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

CREATE TABLE rnaproduct_attrib (

  rnaproduct_id               INT(10) UNSIGNED NOT NULL DEFAULT '0',
  attrib_type_id              SMALLINT(5) UNSIGNED NOT NULL DEFAULT '0',
  value                       TEXT NOT NULL,

  KEY type_val_idx (attrib_type_id, value(40)),
  KEY val_only_idx (value(40)),
  KEY rnaproduct_idx (rnaproduct_id),
  UNIQUE KEY rnaproduct_attribx (rnaproduct_id, attrib_type_id, value(500))

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

CREATE TABLE rnaproduct_type (

  rnaproduct_type_id          SMALLINT(5) UNSIGNED NOT NULL AUTO_INCREMENT,
  code                        VARCHAR(20) NOT NULL DEFAULT '',
  name                        VARCHAR(255) NOT NULL DEFAULT '',
  description                 TEXT,

  PRIMARY KEY (rnaproduct_type_id),
  UNIQUE KEY code_idx (code)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_96_97_c.sql|rnaproduct_tables');
