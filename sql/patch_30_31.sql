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

# gene table changes

ALTER TABLE gene CHANGE type biotype VARCHAR(40) NOT NULL default 'protein_coding';
ALTER TABLE gene ADD source VARCHAR(20) NOT NULL default 'ensembl';
ALTER TABLE gene ADD confidence ENUM( 'KNOWN', 'NOVEL', 'PUTATIVE', 'PREDICTED' ) default 'NOVEL';
ALTER TABLE gene ADD description text;

UPDATE gene g, gene_description gd SET g.description = gd.description WHERE gd.gene_id = g.gene_id;

DROP TABLE gene_description;

# transcript related changes

ALTER TABLE transcript ADD biotype VARCHAR(40) NOT NULL DEFAULT 'protein_coding';
ALTER TABLE transcript ADD confidence ENUM( 'KNOWN', 'NOVEL', 'PUTATIVE', 'PREDICTED' ) default 'NOVEL';
ALTER TABLE transcript ADD description text;

# reasonable default for transcript description
# it might be questionable whether a separate transcript description is necessary
 
UPDATE transcript t, xref x SET t.description = x.description WHERE t.display_xref_id = x.xref_id;
UPDATE transcript SET description=NULL WHERE description="";

# usefull settings for the new tables

# usefull settings for the new tables
UPDATE gene SET biotype='protein_coding' WHERE biotype = 'ensembl';

UPDATE gene g, xref x, external_db ed SET g.confidence='KNOWN' WHERE g.display_xref_id = x.xref_id and x.external_db_id = ed.external_db_id and g.display_xref_id != 0 and ed.status like 'KNOWN%';
UPDATE transcript t, xref x, external_db ed SET t.confidence='KNOWN' WHERE t.display_xref_id = x.xref_id and x.external_db_id = ed.external_db_id and t.display_xref_id != 0 and ed.status like 'KNOWN%';

# some vega specific stuff, shouldnt harm anybody else

UPDATE gene SET biotype='unclassified' WHERE biotype = 'Transcript';
UPDATE gene SET biotype='pseudogene' WHERE biotype = 'Pseudogene';
UPDATE gene SET biotype='protein_coding', confidence='NOVEL' WHERE biotype = 'Novel_CDS';
UPDATE gene SET biotype='unclassified', confidence='NOVEL' WHERE biotype = 'Novel_Transcript';
UPDATE gene SET biotype='unclassified',confidence='PUTATIVE' WHERE biotype = 'Putative';
UPDATE gene SET biotype='protein_coding', confidence='KNOWN' WHERE biotype = 'Known';

UPDATE gene SET biotype='processed_pseudogene' WHERE biotype = 'Processed_pseudogene';
UPDATE gene SET biotype='unprocessed_pseudogene' WHERE biotype = 'Unprocessed_pseudogene';
UPDATE gene SET biotype='protein_coding',confidence='PREDICTED' WHERE biotype = 'Predicted_Gene';
UPDATE gene SET biotype='Ig_segment' WHERE biotype = 'Ig_Segment';
UPDATE gene SET biotype='Ig_pseudogene_segment' WHERE biotype = 'Ig_Pseudogene_Segment';

UPDATE gene SET biotype=replace( biotype, '-','_' );

# reasonable biotypes for the transcripts, take the one from the gene

UPDATE transcript t, gene g SET t.biotype = g.biotype WHERE g.gene_id = t.gene_id;


# new tables regulatory stuff and transcript supporting features

################################################################################
#
# Table structure for table 'regulatory_feature'
#

CREATE TABLE regulatory_feature (

  regulatory_feature_id INT NOT NULL auto_increment,
  name                  VARCHAR(255) NOT NULL,
  seq_region_id         INT NOT NULL,                  # FK refs seq_region
  seq_region_start      INT NOT NULL,
  seq_region_end        INT NOT NULL,
  seq_region_strand     TINYINT NOT NULL,
  analysis_id           INT NOT NULL,                  # FK refs analysis
  regulatory_motif_id   INT,                           # FK refs regulatory_motif
  influence             ENUM('positive', 'negative', 'mixed', 'unknown'),

  PRIMARY KEY(regulatory_feature_id)

) TYPE=MyISAM;

################################################################################
#
# Table structure for table 'regulatory_motif'
#

CREATE TABLE regulatory_motif (

  regulatory_motif_id   INT NOT NULL auto_increment,
  name                  VARCHAR(255) NOT NULL,
  type                  ENUM('miRNA_target', 'promoter'),

  PRIMARY KEY(regulatory_motif_id)

) TYPE=MyISAM;

################################################################################
#
# Table structure for table 'regulatory_feature_object'
#
# Relates regulatory regions to the Ensembl objects they influence. Many-many.

CREATE TABLE regulatory_feature_object (

  regulatory_feature_id INT NOT NULL,               # FK to regulatory_feature
  ensembl_object_type   ENUM( 'Transcript', 'Translation', 'Gene') NOT NULL,
  ensembl_object_id     INT NOT NULL,               # FK to transcript,gene etc

  KEY regulatory_feature_idx (regulatory_feature_id),
  KEY ensembl_object_idx (ensembl_object_type, ensembl_object_id)

) TYPE=MyISAM;


################################################################################
#
# Table structure for table 'peptide_regulatory_feature'
#

CREATE TABLE peptide_regulatory_feature (

  translation_id        INT NOT NULL,               # FK to translation
  regulatory_feature_id INT NOT NULL,               # FK to regulatory_feature

  KEY translation_idx (translation_id),
  KEY regulatory_feature_idx (regulatory_feature_id)

) TYPE=MyISAM;

################################################################################
#
# Table structure for table 'transcript_supporting_feature'
#

CREATE TABLE transcript_supporting_feature (

  transcript_id 	      int(11) DEFAULT '0' NOT NULL,
  feature_type                enum('dna_align_feature','protein_align_feature'),
  feature_id                  int(11) DEFAULT '0' NOT NULL,

  UNIQUE all_idx (transcript_id,feature_type,feature_id),
  KEY feature_idx (feature_type,feature_id)

) MAX_ROWS=100000000 AVG_ROW_LENGTH=80;

