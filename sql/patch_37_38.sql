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

# Patch to convert release 37 Ensembl schema to release 38

UPDATE meta set meta_value="38" where meta_key="schema_version";

# Add info_type & info_text columns to xref

ALTER TABLE xref ADD COLUMN info_type ENUM('PROJECTION', 'MISC');
ALTER TABLE xref ADD COLUMN info_text VARCHAR(255);
ALTER TABLE xref ADD INDEX info_type_idx (info_type);

# Change name of release column in external_db since release is a reserved word
ALTER TABLE external_db CHANGE COLUMN release db_release VARCHAR(40) NOT NULL;

# Add the two new Unmapped Object tables:-


################################################################################
#
# Table structure for table 'unmapped_object'
#
# Describes why a particular external entity was not mapped to an ensembl one.

CREATE TABLE unmapped_object (

  unmapped_object_id    INT UNSIGNED NOT NULL AUTO_INCREMENT,
  type                  ENUM('xref', 'cDNA', 'Marker') NOT NULL,
  analysis_id           INT(10) UNSIGNED NOT NULL,
  external_db_id        INT NOT NULL,
  identifier            VARCHAR(255) NOT NULL,
  unmapped_reason_id    SMALLINT(5) UNSIGNED NOT NULL,
  query_score           DOUBLE,
  target_score          DOUBLE,
  ensembl_id            INT(10) unsigned default '0',
  ensembl_object_type   ENUM('RawContig','Transcript','Gene','Translation') collate latin1_bin default 'RawContig',
  PRIMARY KEY            ( unmapped_object_id ),
  KEY                    id_idx( identifier ),
  KEY                    anal_idx( analysis_id ),
  KEY                    anal_exdb_idx( analysis_id, external_db_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;

################################################################################
#
# Table structure for table 'unmapped_reason'
#
# Describes the reason why a mapping failed.

CREATE TABLE unmapped_reason (

  unmapped_reason_id     SMALLINT(5) UNSIGNED NOT NULL AUTO_INCREMENT,
  summary_description    VARCHAR(255),
  full_description       VARCHAR(255),

  PRIMARY KEY ( unmapped_reason_id )

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


# Add some more columns to mapping_session

ALTER TABLE mapping_session ADD COLUMN new_assembly varchar(20) NOT NULL default '' AFTER new_db_name;
ALTER TABLE mapping_session ADD COLUMN old_assembly varchar(20) NOT NULL default '' AFTER new_db_name;
ALTER TABLE mapping_session ADD COLUMN new_release varchar(5) NOT NULL default '' AFTER new_db_name;
ALTER TABLE mapping_session ADD COLUMN old_release varchar(5) NOT NULL default '' AFTER new_db_name;
ALTER TABLE mapping_session CHANGE created created DATETIME NOT NULL;


# Add the new oligo tables, which replace the affy tables

CREATE TABLE oligo_feature (
       oligo_feature_id   INT NOT NULL auto_increment,
       seq_region_id      INT UNSIGNED NOT NULL,
       seq_region_start   INT NOT NULL,
       seq_region_end     INT NOT NULL,
       seq_region_strand  TINYINT NOT NULL,
       
       mismatches         TINYINT,
       oligo_probe_id     INT NOT NULL,
       analysis_id        INT(10) UNSIGNED NOT NULL,

       PRIMARY KEY (oligo_feature_id),
       KEY seq_region_idx (seq_region_id, seq_region_start),
       KEY probe_idx (oligo_probe_id)
) COLLATE=latin1_swedish_ci TYPE=MyISAM;

CREATE TABLE oligo_probe (
       oligo_probe_id  INT NOT NULL auto_increment,
       oligo_array_id  INT NOT NULL,
       probeset        VARCHAR(40),
       name            VARCHAR(20),
       description     TEXT,
       length          SMALLINT NOT NULL,

       PRIMARY KEY (oligo_probe_id, oligo_array_id),
       KEY probeset_idx (probeset),
       KEY array_idx (oligo_array_id)
) COLLATE=latin1_swedish_ci TYPE=MyISAM;

CREATE TABLE oligo_array (
       oligo_array_id   INT NOT NULL auto_increment,
       parent_array_id  INT,
       probe_setsize    TINYINT NOT NULL,
       name             VARCHAR(40) NOT NULL,
       type             ENUM( 'AFFY', 'OLIGO' ),

       PRIMARY KEY (oligo_array_id)
) COLLATE=latin1_swedish_ci TYPE=MyISAM;


# Copy any existing data from the affy tables to the oligo tables

INSERT INTO oligo_array (oligo_array_id, parent_array_id,
       probe_setsize, name, type)
SELECT affy_array_id, parent_array_id, probe_setsize, name, 'AFFY'
FROM   affy_array;

INSERT INTO oligo_probe (oligo_probe_id, oligo_array_id,
       probeset, name, description, length)
SELECT affy_probe_id, affy_array_id, probeset, name, NULL, 25
FROM   affy_probe;

INSERT INTO oligo_feature (oligo_feature_id, seq_region_id,
       seq_region_start, seq_region_end, seq_region_strand,
       mismatches, oligo_probe_id, analysis_id)
SELECT affy_feature_id, seq_region_id, seq_region_start, seq_region_end,
       seq_region_strand, mismatches, affy_probe_id, analysis_id
FROM   affy_feature;

UPDATE meta_coord
SET    table_name='oligo_feature'
WHERE  table_name='affy_feature';

# Drop the affy tables

DROP TABLE affy_array;
DROP TABLE affy_probe;
DROP TABLE affy_feature;


# change column 'value' in attrib tables from VARCHAR(255) to TEXT

ALTER TABLE misc_attrib DROP INDEX type_val_idx;
ALTER TABLE misc_attrib CHANGE value value TEXT NOT NULL DEFAULT '';
ALTER TABLE misc_attrib ADD INDEX type_val_idx (attrib_type_id, value(40));

ALTER TABLE seq_region_attrib DROP INDEX type_val_idx;
ALTER TABLE seq_region_attrib CHANGE value value TEXT NOT NULL DEFAULT '';
ALTER TABLE seq_region_attrib ADD INDEX type_val_idx (attrib_type_id, value(40));

ALTER TABLE gene_attrib DROP INDEX type_val_idx;
ALTER TABLE gene_attrib CHANGE value value TEXT NOT NULL DEFAULT '';
ALTER TABLE gene_attrib ADD INDEX type_val_idx (attrib_type_id, value(40));

ALTER TABLE transcript_attrib DROP INDEX type_val_idx;
ALTER TABLE transcript_attrib CHANGE value value TEXT NOT NULL DEFAULT '';
ALTER TABLE transcript_attrib ADD INDEX type_val_idx (attrib_type_id, value(40));

ALTER TABLE translation_attrib DROP INDEX type_val_idx;
ALTER TABLE translation_attrib CHANGE value value TEXT NOT NULL DEFAULT '';
ALTER TABLE translation_attrib ADD INDEX type_val_idx (attrib_type_id, value(40));


# consistent analysis_id foreign key columns across database

ALTER TABLE gene CHANGE analysis_id analysis_id INT(10) UNSIGNED NOT NULL;
ALTER TABLE identity_xref CHANGE analysis_id analysis_id INT(10) UNSIGNED NOT NULL;
ALTER TABLE prediction_transcript CHANGE analysis_id analysis_id INT(10) UNSIGNED NOT NULL;
ALTER TABLE qtl_feature CHANGE analysis_id analysis_id INT(10) UNSIGNED NOT NULL;
ALTER TABLE density_type CHANGE analysis_id analysis_id INT(10) UNSIGNED NOT NULL;
ALTER TABLE regulatory_feature CHANGE analysis_id analysis_id INT(10) UNSIGNED NOT NULL;
ALTER TABLE regulatory_search_region CHANGE analysis_id analysis_id INT(10) UNSIGNED NOT NULL;


# add column analysis_id to transcript

ALTER TABLE transcript ADD COLUMN analysis_id INT(10) UNSIGNED NOT NULL AFTER gene_id;
ALTER TABLE transcript ADD INDEX analysis_idx (analysis_id);
UPDATE transcript t, gene g set t.analysis_id = g.analysis_id WHERE t.gene_id = g.gene_id;

