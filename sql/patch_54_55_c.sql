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

# patch_54_55_c.sql
#
# title: Add alternative splicing event tables
#
# description:
# Add alternative splicing event tables for annotation

CREATE TABLE splicing_event (

  splicing_event_id       INT(10)  UNSIGNED NOT NULL AUTO_INCREMENT,
  name                    VARCHAR(134),
  gene_id                 INT(10) UNSIGNED NOT NULL,
  seq_region_id           INT(10) UNSIGNED NOT NULL,
  seq_region_start        INT(10) UNSIGNED NOT NULL,
  seq_region_end          INT(10) UNSIGNED NOT NULL,
  seq_region_strand       TINYINT(2) NOT NULL,
  type	                  ENUM('CNE','CE','AFE','A5SS','A3SS','MXE','IR','II','EI', 'AT', 'ALE', 'AI'),
  PRIMARY KEY (splicing_event_id),
  KEY gene_idx (gene_id),
  KEY seq_region_idx (seq_region_id, seq_region_start)

)  COLLATE=latin1_swedish_ci TYPE=MyISAM;


CREATE TABLE splicing_event_feature (

  splicing_event_feature_id INT(10)  UNSIGNED NOT NULL,
  splicing_event_id         INT(10)  UNSIGNED NOT NULL,
  exon_id                   INT(10)  UNSIGNED NOT NULL,
  transcript_id             INT(10)  UNSIGNED NOT NULL,
  feature_order             INT(10)  UNSIGNED NOT NULL,
  transcript_association    INT(10)  UNSIGNED NOT NULL,
  type                      ENUM('constitutive_exon','exon','flanking_exon'),
  start                     INT(10)  UNSIGNED NOT NULL,
  end                       INT(10)  UNSIGNED NOT NULL,

  PRIMARY KEY (splicing_event_feature_id,exon_id,transcript_id),
  KEY se_idx (splicing_event_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;



CREATE TABLE splicing_transcript_pair (


  splicing_transcript_pair_id INT(10)  UNSIGNED NOT NULL,
  splicing_event_id           INT(10)  UNSIGNED NOT NULL, 
  transcript_id_1             INT(10)  UNSIGNED NOT NULL,
  transcript_id_2             INT(10)  UNSIGNED NOT NULL,

  PRIMARY KEY (splicing_transcript_pair_id),
  KEY se_idx (splicing_event_id)
  

) COLLATE=latin1_swedish_ci TYPE=MyISAM;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_54_55_c.sql|add_splicing_event_tables');
