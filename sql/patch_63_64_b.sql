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

