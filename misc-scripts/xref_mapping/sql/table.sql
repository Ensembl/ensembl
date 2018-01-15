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

# Schema for internal-external database mappings (xrefs)


################################################################################
#
# General external annotation.

CREATE TABLE xref (

  xref_id                     int unsigned not null auto_increment,
  accession                   varchar(255) not null,
  version                     int unsigned,
  label                       varchar(255),
  description                 text,
  source_id                   int unsigned not null,
  species_id                  int unsigned not null,
  info_type                   ENUM( 'NONE', 'PROJECTION', 'MISC', 'DEPENDENT',
                                    'DIRECT', 'SEQUENCE_MATCH',
                                    'INFERRED_PAIR', 'PROBE',
                                    'UNMAPPED', 'COORDINATE_OVERLAP',
                                    'CHECKSUM') DEFAULT 'NONE' NOT NULL,
  info_text	              VARCHAR(255) DEFAULT '' NOT NULL,
  dumped                      ENUM( 'MAPPED', 'NO_DUMP_ANOTHER_PRIORITY', 'UNMAPPED_NO_MAPPING',
                                    'UNMAPPED_NO_MASTER', 'UNMAPPED_MASTER_FAILED', 
                                    'UNMAPPED_NO_STABLE_ID', 'UNMAPPED_INTERPRO') DEFAULT null,

  PRIMARY KEY (xref_id),
  UNIQUE acession_idx(accession,label,source_id,species_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

################################################################################

CREATE TABLE primary_xref (

  xref_id                     int unsigned not null,
  sequence                    mediumtext,
  sequence_type               enum('dna','peptide'),
  status                      enum('experimental','predicted'),

  PRIMARY KEY (xref_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

################################################################################

CREATE TABLE dependent_xref (

  object_xref_id              int unsigned,
  master_xref_id              int unsigned not null,
  dependent_xref_id           int unsigned not null,
  linkage_annotation          varchar(255),
  linkage_source_id           int unsigned not null,

  KEY master_idx(master_xref_id),
  KEY dependent_idx(dependent_xref_id),
  KEY object_id(object_xref_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;





################################################################################

CREATE TABLE synonym (

  xref_id                     int unsigned not null,
  synonym                     varchar(255),

  UNIQUE KEY (xref_id, synonym),
  KEY xref_idx(xref_id),
  KEY synonym_idx(synonym)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

################################################################################

CREATE TABLE dependent_source (
  master_source_id           int unsigned not null,
  dependent_name             varchar(255) not null,

  PRIMARY KEY (master_source_id, dependent_name)
) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

################################################################################

CREATE TABLE source (

  source_id                   int unsigned not null auto_increment,
  name                        varchar(255) not null,
  status                      enum('KNOWN','XREF','PRED','ORTH','PSEUDO','LOWEVIDENCE','NOIDEA') not null default 'NOIDEA',
  source_release              varchar(255),
  download                    enum('Y','N') default 'Y',
  ordered                     int unsigned not null, 
  priority                    int unsigned default 1,
  priority_description        varchar(40) default "",
   
  PRIMARY KEY (source_id),
  KEY name_idx(name) 

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

################################################################################

CREATE TABLE source_url (

  source_url_id               int unsigned not null auto_increment,
  source_id                   int unsigned not null,
  species_id                  int unsigned not null,
  url                         mediumtext,
  release_url                 mediumtext,
  checksum                    varchar(1025),
  file_modified_date          datetime,
  upload_date                 datetime,
  parser                      varchar(255),

  PRIMARY KEY (source_url_id),
  KEY source_idx(source_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

################################################################################


CREATE TABLE source_mapping_method (
       
       source_id                   int unsigned not null,
       method         VARCHAR(255),   

       UNIQUE KEY (source_id, method)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


CREATE TABLE gene_direct_xref (

  general_xref_id             int unsigned not null,
  ensembl_stable_id           varchar(255),
  linkage_xref                varchar(255),

  KEY primary_idx(general_xref_id),
  KEY ensembl_idx(ensembl_stable_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


CREATE TABLE transcript_direct_xref (

  general_xref_id             int unsigned not null,
  ensembl_stable_id           varchar(255),
  linkage_xref                varchar(255),

  KEY primary_idx(general_xref_id),
  KEY ensembl_idx(ensembl_stable_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

CREATE TABLE translation_direct_xref (

  general_xref_id             int unsigned not null,
  ensembl_stable_id           varchar(255),
  linkage_xref                varchar(255),

  KEY primary_idx(general_xref_id),
  KEY ensembl_idx(ensembl_stable_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

################################################################################

CREATE TABLE species (

  species_id                  int unsigned not null,
  taxonomy_id                 int unsigned not null,
  name                        varchar(255) not null,
  aliases                     varchar(255),

  KEY species_idx (species_id),
  KEY taxonomy_idx(taxonomy_id),
  UNIQUE KEY species_taxonomy_idx(species_id,taxonomy_id),
  KEY name_idx(name)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

################################################################################

CREATE TABLE interpro (

  interpro               varchar(255) not null,
  pfam                   varchar(255) not null,
  dbtype                 enum ('PROSITE','PFAM','PREFILE','PROFILE','TIGRFAMs','PRINTS','PIRSF','SMART','SSF')  not null

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

################################################################################

CREATE TABLE pairs (

  source_id			 int unsigned not null,
  accession1                     varchar(255) not null,
  accession2                     varchar(255) not null,

  KEY ac2_idx(accession2)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;
################################################################################

-- Table for coordinate-based Xrefs, based
-- on the knownGenes table from UCSC.

CREATE TABLE coordinate_xref (
  coord_xref_id INT UNSIGNED NOT NULL AUTO_INCREMENT,
  source_id     INT UNSIGNED NOT NULL,
  species_id    INT UNSIGNED NOT NULL,
  accession     VARCHAR(255) NOT NULL,
  chromosome    VARCHAR(255) NOT NULL,
  strand        TINYINT(2) NOT NULL,
  txStart       INT(10) NOT NULL,
  txEnd         INT(10) NOT NULL,
  cdsStart      INT(10),
  cdsEnd        INT(10),
  exonStarts    TEXT NOT NULL,
  exonEnds      TEXT NOT NULL,

  UNIQUE KEY coord_xref_idx(coord_xref_id),
  INDEX start_pos_idx(species_id, chromosome, strand, txStart),
  INDEX end_pos_idx(species_id, chromosome, strand, txEnd)
) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

################################################################################

-- Table for checksum-based Xrefs, based
-- on the input format from UniProt/UniParc. This is MyISAM because
-- we do a LOAD DATA INFILE & the MyISAM engine does not handle this too
-- well (especially since there is no need for transactions for this table)

CREATE TABLE checksum_xref (
  checksum_xref_id  INT UNSIGNED NOT NULL AUTO_INCREMENT,
  source_id         INT UNSIGNED NOT NULL,
  accession         CHAR(14) NOT NULL,
  checksum          CHAR(32) NOT NULL,

  PRIMARY KEY (checksum_xref_id),
  INDEX checksum_idx(checksum(10))
) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


################################################################################

################################################################################

-- new tables for new mapper code

CREATE TABLE mapping (
  job_id         INT UNSIGNED,
  type           enum('dna','peptide','UCSC'), # not sure about UCSC
  command_line   text,
  percent_query_cutoff    INT UNSIGNED,
  percent_target_cutoff   INT UNSIGNED,
  method         VARCHAR(255),
  array_size     INT UNSIGNED

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

CREATE TABLE mapping_jobs (
  root_dir          text,
  map_file          VARCHAR(255),
  status            enum('SUBMITTED','FAILED','SUCCESS'),
  out_file          VARCHAR(255),
  err_file          VARCHAR(255),
  array_number      INT UNSIGNED,
  job_id            INT UNSIGNED,
  failed_reason     VARCHAR(255),
  object_xref_start INT UNSIGNED,
  object_xref_end   INT UNSIGNED

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

CREATE TABLE gene_transcript_translation (

  gene_id			INT UNSIGNED NOT NULL,
  transcript_id			INT UNSIGNED NOT NULL,
  translation_id		INT UNSIGNED,
  PRIMARY KEY (transcript_id),
  INDEX gene_idx (gene_id),
  INDEX translation_idx (translation_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

CREATE TABLE core_database (
  port		INT UNSIGNED,
  user          VARCHAR(16),
  pass          VARCHAR(16),
  dbname        VARCHAR(16),
  xref_dir      text,
  core_dir      text

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;
  


CREATE TABLE havana_status (

  stable_id    VARCHAR(128),
  status       enum('KNOWN','NOVEL','PUTATIVE','PREDICTED','KNOWN_BY_PROJECTION','UNKNOWN'),
  UNIQUE KEY status_idx(stable_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


-- Try to keep the status in the correct order
--   it will make it easier to see what is happening


CREATE TABLE process_status (
  id            INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  status	enum('xref_created','parsing_started','parsing_finished','alt_alleles_added',
                     'xref_fasta_dumped','core_fasta_dumped','core_data_loaded',
                     'mapping_submitted','mapping_finished','mapping_processed',
                     'direct_xrefs_parsed',
                     'prioritys_flagged','processed_pairs','biomart_test_finished',
                     'source_level_move_finished','alt_alleles_processed',	
                     'official_naming_done',
                     'checksum_xrefs_started', 'checksum_xrefs_finished',
                     'coordinate_xrefs_started','coordinate_xref_finished',
                     'tests_started','tests_failed','tests_finished',
                     'core_loaded','display_xref_done','gene_description_done'),
  date          DATETIME NOT NULL,
  PRIMARY KEY (id)
) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


-- This table is populated in DisplayXrefs.pm and is used to set the best xrefs as 
-- display xrefs for genes and transcripts



CREATE TABLE display_xref_priority(
    ensembl_object_type	VARCHAR(100),
    source_id INT UNSIGNED NOT NULL,
    priority  SMALLINT UNSIGNED NOT NULL,

    PRIMARY KEY (ensembl_object_type, source_id)
) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


-- This table is populated in DisplayXrefs.pm and is used to set
-- gene descriptions

CREATE TABLE gene_desc_priority(
    source_id	   INT UNSIGNED NOT NULL,
    priority       SMALLINT UNSIGNED NOT NULL,

    PRIMARY KEY (source_id)
) COLLATE=latin1_swedish_ci ENGINE=MyISAM;
                     

################################################################################
################################################################################
################################################################################

-- Incorporated but modified core tables

CREATE TABLE alt_allele (
  alt_allele_id         INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  gene_id               INT(10) UNSIGNED NOT NULL,
  is_reference          INT UNSIGNED DEFAULT 0,

  UNIQUE KEY gene_idx (gene_id),
  UNIQUE KEY allele_idx (alt_allele_id, gene_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;



CREATE TABLE gene_stable_id (

  internal_id                 INT UNSIGNED NOT NULL,
  stable_id                   VARCHAR(128) NOT NULL,
  display_xref_id             INT UNSIGNED DEFAULT NULL,
  desc_set                    INT UNSIGNED DEFAULT 0,               

  PRIMARY KEY (stable_id),
  INDEX internal_idx (internal_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

CREATE TABLE transcript_stable_id (

  internal_id                 INT UNSIGNED NOT NULL,
  stable_id                   VARCHAR(128) NOT NULL,
  display_xref_id             INT UNSIGNED DEFAULT NULL,
  biotype                     VARCHAR(40) NOT NULL,

  PRIMARY KEY (stable_id),
  INDEX internal_idx (internal_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

CREATE TABLE translation_stable_id (

  internal_id                 INT UNSIGNED NOT NULL,
  stable_id                   VARCHAR(128) NOT NULL,

  PRIMARY KEY (internal_id),
  INDEX stable_idx (stable_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;




CREATE TABLE object_xref (

  object_xref_id              INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  ensembl_id                  INT(10) UNSIGNED NOT NULL,
  ensembl_object_type         ENUM('RawContig', 'Transcript', 'Gene',
                                   'Translation')
                              NOT NULL,
  xref_id                     INT UNSIGNED NOT NULL,
  linkage_annotation          VARCHAR(255) DEFAULT NULL,
  linkage_type                ENUM( 'PROJECTION', 'MISC', 'DEPENDENT',
                                    'DIRECT', 'SEQUENCE_MATCH',
                                    'INFERRED_PAIR', 'PROBE',
                                    'UNMAPPED', 'COORDINATE_OVERLAP',
                                    'CHECKSUM'),
  ox_status                   ENUM( 'DUMP_OUT','FAILED_PRIORITY', 'FAILED_CUTOFF', 'NO_DISPLAY', 'MULTI_DELETE')  NOT NULL DEFAULT 'DUMP_OUT',
-- set ox_status to 0 if non used priority_xref or failed cutoff
  unused_priority             INT UNSIGNED,
  master_xref_id              INT UNSIGNED DEFAULT NULL,

  PRIMARY KEY (object_xref_id),
  UNIQUE (ensembl_object_type, ensembl_id, xref_id, ox_status, master_xref_id),
  KEY oxref_idx (object_xref_id, xref_id, ensembl_object_type, ensembl_id),
  KEY xref_idx (xref_id, ensembl_object_type)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


CREATE TABLE identity_xref (

  object_xref_id          INT(10) UNSIGNED NOT NULL,
  query_identity          INT(5),
  target_identity         INT(5),

  hit_start               INT,
  hit_end                 INT,
  translation_start       INT,
  translation_end         INT,
  cigar_line              TEXT,

  score                   DOUBLE,
  evalue                  DOUBLE,
--  analysis_id             SMALLINT UNSIGNED NOT NULL, # set in core not needed in xref

  PRIMARY KEY (object_xref_id)
--  KEY analysis_idx (analysis_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

CREATE TABLE go_xref (

  object_xref_id          INT(10) UNSIGNED DEFAULT '0' NOT NULL,
  linkage_type            CHAR(3) NOT NULL,
  source_xref_id          INT(10) UNSIGNED DEFAULT NULL,
  KEY (object_xref_id),
  KEY (source_xref_id),
  UNIQUE (object_xref_id, source_xref_id, linkage_type)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


CREATE TABLE meta (

  meta_id                     INT NOT NULL AUTO_INCREMENT,
  species_id                  INT UNSIGNED DEFAULT 1,
  meta_key                    VARCHAR(40) NOT NULL,
  meta_value                  VARCHAR(255) BINARY NOT NULL,
  date                        DATETIME NOT NULL,

  PRIMARY   KEY (meta_id),
  UNIQUE    KEY species_key_value_idx (meta_id, species_id, meta_key, meta_value),
            KEY species_value_idx (species_id, meta_value)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;




