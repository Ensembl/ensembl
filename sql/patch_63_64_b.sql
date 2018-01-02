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

# patch_63_64_b.sql
#
# Title: Add tables for supporting operons
#
# Description:
# Create operon, operon_transcript, operon_transcript_gene, operon_stable_id, operon_transcript_stable_id tables to represent Operon and OperonTranscript and to associated them with genes. Modify object_xref to allow xrefs to be attached to operon and operon_transcript

CREATE TABLE operon (

  operon_id                     INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  seq_region_id               INT(10) UNSIGNED NOT NULL,
  seq_region_start            INT(10) UNSIGNED NOT NULL,
  seq_region_end              INT(10) UNSIGNED NOT NULL,
  seq_region_strand           TINYINT(2) NOT NULL,
  display_label        		      VARCHAR(255) DEFAULT NULL,
  analysis_id                 SMALLINT UNSIGNED NOT NULL,
  PRIMARY KEY (operon_id),
  KEY seq_region_idx (seq_region_id, seq_region_start),
  KEY name_idx (display_label)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

CREATE TABLE operon_transcript (

  operon_transcript_id        INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  seq_region_id               INT(10) UNSIGNED NOT NULL,
  seq_region_start            INT(10) UNSIGNED NOT NULL,
  seq_region_end              INT(10) UNSIGNED NOT NULL,
  seq_region_strand           TINYINT(2) NOT NULL,
  operon_id                     INT(10) UNSIGNED NOT NULL,
  display_label        		      VARCHAR(255) DEFAULT NULL,
  analysis_id                 SMALLINT UNSIGNED NOT NULL,
  PRIMARY KEY (operon_transcript_id),
  KEY operon_idx (operon_id),
  KEY seq_region_idx (seq_region_id, seq_region_start)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

CREATE TABLE operon_transcript_gene (

  operon_transcript_id        INT(10) UNSIGNED,
  gene_id                     INT(10) UNSIGNED,

  KEY operon_transcript_gene_idx (operon_transcript_id,gene_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

CREATE TABLE operon_stable_id (

  operon_id                     INT UNSIGNED NOT NULL,
  stable_id                   VARCHAR(128) NOT NULL,
  version                     INT(10),
  created_date                DATETIME NOT NULL,
  modified_date               DATETIME NOT NULL,

  PRIMARY KEY (operon_id),
  KEY stable_id_idx (stable_id, version)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

CREATE TABLE operon_transcript_stable_id (

  operon_transcript_id                     INT UNSIGNED NOT NULL,
  stable_id                   VARCHAR(128) NOT NULL,
  version                     INT(10),
  created_date                DATETIME NOT NULL,
  modified_date               DATETIME NOT NULL,

  PRIMARY KEY (operon_transcript_id),
  KEY stable_id_idx (stable_id, version)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

alter table object_xref modify column ensembl_object_type ENUM('RawContig', 'Transcript', 'Gene', 'Translation', 'Operon', 'OperonTranscript') NOT NULL;

# Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_63_64_b.sql|add_operons');

