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

# patch_92_93_b.sql
#
# Title: Added biotype table
#
# Description:
#   Added new table biotype

CREATE TABLE biotype (
  biotype_id      INTEGER UNSIGNED NOT NULL AUTO_INCREMENT,
  name            VARCHAR(64) NOT NULL,
  object_type     ENUM('gene','transcript') NOT NULL DEFAULT 'gene',
  db_type         set('cdna','core','coreexpressionatlas','coreexpressionest','coreexpressiongnf','funcgen','otherfeatures','rnaseq','variation','vega','presite','sangervega') NOT NULL DEFAULT 'core',
  attrib_type_id  INTEGER DEFAULT NULL,
  description     TEXT,
  biotype_group   ENUM('coding','pseudogene','snoncoding','lnoncoding','mnoncoding','LRG','undefined','no_group') DEFAULT NULL,
  so_acc          VARCHAR(64),
  PRIMARY KEY (biotype_id),
  UNIQUE KEY name_type_idx (name, object_type)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_92_93_b.sql|biotype_table');


