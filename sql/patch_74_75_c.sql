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

# patch_74_75_c.sql
#
# title: Add genome_statistics table
#
# description:
# Addition of a new table, genome_statistics, to store genome related statistics

CREATE TABLE genome_statistics (

  genome_statistics_id     INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  statistic                VARCHAR(128) NOT NULL,
  value                    INT(10) UNSIGNED DEFAULT '0' NOT NULL,
  species_id               INT UNSIGNED DEFAULT 1,
  attrib_type_id           INT(10) UNSIGNED  DEFAULT NULL,
  timestamp                DATETIME NOT NULL DEFAULT '0000-00-00 00:00:00',

  PRIMARY KEY (genome_statistics_id),
  UNIQUE KEY stats_uniq(statistic, attrib_type_id, species_id),
  KEY stats_idx (statistic, attrib_type_id, species_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_74_75_c.sql|add_genome_statistics');

 
