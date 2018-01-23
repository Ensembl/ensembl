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

# patch_64_65_b.sql
#
# Title: Incorporate stable_id tables into exon, gene, operon, operon_transcript, translation and transcript tables.
#
# Description: 

# Add stable_id, version, stable_id_created_date, stable_id_modified_date columns to tables exon, 
# gene, operon, operon_transcript, transcript and translation. Add index on stable_id and version to each table.
# Move data from the stable_id tables into their related object tables.
# Delete stable_id tables. Create views like the stable_id tables.


ALTER TABLE exon ADD COLUMN (stable_id VARCHAR(128) DEFAULT NULL, version SMALLINT UNSIGNED NOT NULL DEFAULT 1, created_date DATETIME NOT NULL DEFAULT '0000-00-00 00:00:00', modified_date DATETIME NOT NULL DEFAULT '0000-00-00 00:00:00');
ALTER TABLE exon ADD INDEX stable_id_idx (stable_id, version);

ALTER TABLE gene ADD COLUMN (stable_id VARCHAR(128) DEFAULT NULL, version SMALLINT UNSIGNED NOT NULL DEFAULT 1, created_date DATETIME NOT NULL DEFAULT '0000-00-00 00:00:00', modified_date DATETIME NOT NULL DEFAULT '0000-00-00 00:00:00');
ALTER TABLE gene ADD INDEX stable_id_idx (stable_id, version);

ALTER TABLE operon ADD COLUMN (stable_id VARCHAR(128) DEFAULT NULL, version SMALLINT UNSIGNED NOT NULL DEFAULT 1, created_date DATETIME NOT NULL DEFAULT '0000-00-00 00:00:00', modified_date DATETIME NOT NULL DEFAULT '0000-00-00 00:00:00');
ALTER TABLE operon ADD INDEX stable_id_idx (stable_id, version);

ALTER TABLE operon_transcript ADD COLUMN (stable_id VARCHAR(128) DEFAULT NULL, version SMALLINT UNSIGNED NOT NULL DEFAULT 1, created_date DATETIME NOT NULL DEFAULT '0000-00-00 00:00:00', modified_date DATETIME NOT NULL DEFAULT '0000-00-00 00:00:00');
ALTER TABLE operon_transcript ADD INDEX stable_id_idx (stable_id, version);

ALTER TABLE transcript ADD COLUMN (stable_id VARCHAR(128) DEFAULT NULL, version SMALLINT UNSIGNED NOT NULL DEFAULT 1, created_date DATETIME NOT NULL DEFAULT '0000-00-00 00:00:00', modified_date DATETIME NOT NULL DEFAULT '0000-00-00 00:00:00');
ALTER TABLE transcript ADD INDEX stable_id_idx (stable_id, version);
  
ALTER TABLE translation ADD COLUMN (stable_id VARCHAR(128) DEFAULT NULL, version SMALLINT UNSIGNED NOT NULL DEFAULT 1, created_date DATETIME NOT NULL DEFAULT '0000-00-00 00:00:00', modified_date DATETIME NOT NULL DEFAULT '0000-00-00 00:00:00');
ALTER TABLE translation ADD INDEX stable_id_idx (stable_id, version);

UPDATE exon o, exon_stable_id s SET o.stable_id = s.stable_id, o.version = s.version, o.created_date = s.created_date, o.modified_date = s.modified_date WHERE o.exon_id = s.exon_id;

UPDATE gene o, gene_stable_id s SET o.stable_id = s.stable_id, o.version = s.version, o.created_date = s.created_date, o.modified_date = s.modified_date WHERE o.gene_id = s.gene_id;

UPDATE operon o, operon_stable_id s SET o.stable_id = s.stable_id, o.version = s.version, o.created_date = s.created_date, o.modified_date = s.modified_date WHERE o.operon_id = s.operon_id;

UPDATE operon_transcript o, operon_transcript_stable_id s SET o.stable_id = s.stable_id, o.version = s.version, o.created_date = s.created_date, o.modified_date = s.modified_date WHERE o.operon_transcript_id = s.operon_transcript_id;

UPDATE translation o, translation_stable_id s SET o.stable_id = s.stable_id, o.version = s.version, o.created_date = s.created_date, o.modified_date = s.modified_date WHERE o.translation_id = s.translation_id;

UPDATE transcript o, transcript_stable_id s SET o.stable_id = s.stable_id, o.version = s.version, o.created_date = s.created_date, o.modified_date = s.modified_date WHERE o.transcript_id = s.transcript_id;

DROP TABLE exon_stable_id;

DROP TABLE gene_stable_id;

DROP TABLE operon_stable_id;

DROP TABLE operon_transcript_stable_id;

DROP TABLE translation_stable_id;

DROP TABLE transcript_stable_id;

CREATE DEFINER = CURRENT_USER SQL SECURITY INVOKER VIEW exon_stable_id (exon_id, stable_id, version, created_date, modified_date) AS (SELECT exon_id, stable_id, version, created_date, modified_date FROM exon);

CREATE DEFINER = CURRENT_USER SQL SECURITY INVOKER VIEW gene_stable_id (gene_id, stable_id, version, created_date, modified_date) AS (SELECT gene_id, stable_id, version, created_date, modified_date FROM gene);

CREATE DEFINER = CURRENT_USER SQL SECURITY INVOKER VIEW operon_stable_id (operon_id, stable_id, version, created_date, modified_date) AS (SELECT operon_id, stable_id, version, created_date, modified_date FROM operon);

CREATE DEFINER = CURRENT_USER SQL SECURITY INVOKER VIEW operon_transcript_stable_id (operon_transcript_id, stable_id, version, created_date, modified_date) AS (SELECT operon_transcript_id, stable_id, version, created_date, modified_date FROM operon_transcript);

CREATE DEFINER = CURRENT_USER SQL SECURITY INVOKER VIEW translation_stable_id (translation_id, stable_id, version, created_date, modified_date) AS (SELECT translation_id, stable_id, version, created_date, modified_date FROM translation);

CREATE DEFINER = CURRENT_USER SQL SECURITY INVOKER VIEW transcript_stable_id (transcript_id, stable_id, version, created_date, modified_date) AS (SELECT transcript_id, stable_id, version, created_date, modified_date FROM transcript);


# Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_64_65_b.sql|merge_stable_id_with_object');

