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

# Add linkable & priority columns to external_db
ALTER TABLE external_db ADD COLUMN dbprimary_acc_linkable BOOLEAN DEFAULT 1 NOT NULL;
ALTER TABLE external_db ADD COLUMN display_label_linkable BOOLEAN DEFAULT 0 NOT NULL;
ALTER TABLE external_db ADD COLUMN priority INT NOT NULL;

################################################################################

# Many changes to regulatory-feature-related tables
# Since these tables are empty it is easier to drop the existing 
# tables and re-create the new ones

DROP TABLE regulatory_feature;
DROP TABLE regulatory_motif;
DROP TABLE regulatory_feature_object;
DROP TABLE peptide_regulatory_feature;

################################################################################
#
# Table structure for table 'regulatory_feature'
#
# Describes instances of regulatory_factor binding to the genome.

CREATE TABLE regulatory_feature (

  regulatory_feature_id INT NOT NULL auto_increment,
  name                  VARCHAR(255) NOT NULL,
  seq_region_id         INT NOT NULL,                  # FK refs seq_region
  seq_region_start      INT NOT NULL,
  seq_region_end        INT NOT NULL,
  seq_region_strand     TINYINT NOT NULL,
  analysis_id           INT NOT NULL,                  # FK refs analysis
  regulatory_factor_id  INT,                           # FK refs regulatory_motif


  PRIMARY KEY(regulatory_feature_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;

################################################################################
#
# Table structure for table 'regulatory_factor'
#

CREATE TABLE regulatory_factor (

  regulatory_factor_id   INT NOT NULL auto_increment,
  name                   VARCHAR(255) NOT NULL,
  type                   ENUM('miRNA_target', 'transcription_factor', 'transcription_factor_complex'),

  PRIMARY KEY(regulatory_factor_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;

################################################################################
#
# Table structure for table 'regulatory_feature_object'
#
# Relates regulatory regions to the Ensembl objects they influence. Many-many.

CREATE TABLE regulatory_feature_object (

  regulatory_feature_id INT NOT NULL,               # FK to regulatory_feature
  ensembl_object_type   ENUM( 'Transcript', 'Translation', 'Gene') NOT NULL,
  ensembl_object_id     INT NOT NULL,               # FK to transcript,gene etc
  influence             ENUM('positive', 'negative', 'mixed', 'unknown'),
  evidence              VARCHAR(255),

  KEY regulatory_feature_idx (regulatory_feature_id),
  KEY ensembl_object_idx (ensembl_object_type, ensembl_object_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'regulatory_factor_transcript'
#
# Describes which transcripts code for particular regulatory factors.

CREATE TABLE regulatory_factor_transcript (

  transcript_id         INT NOT NULL,               # FK to transcript
  regulatory_factor_id  INT NOT NULL,               # FK to regulatory_feature

  KEY translation_idx (transcript_id),
  KEY regulatory_factor_idx (regulatory_factor_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;

