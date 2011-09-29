CREATE TABLE transcript (

  transcript_id               INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  gene_id                     INT(10) UNSIGNED,
  analysis_id                 SMALLINT UNSIGNED NOT NULL,
  seq_region_id               INT(10) UNSIGNED NOT NULL,
  seq_region_start            INT(10) UNSIGNED NOT NULL,
  seq_region_end              INT(10) UNSIGNED NOT NULL,
  seq_region_strand           TINYINT(2) NOT NULL,
  display_xref_id             INT(10) UNSIGNED,
  biotype                     VARCHAR(40) NOT NULL,
  status                      ENUM('KNOWN', 'NOVEL', 'PUTATIVE', 'PREDICTED', 'KNOWN_BY_PROJECTION', 'UNKNOWN'),
  description                 TEXT,
  is_current                  BOOLEAN NOT NULL DEFAULT 1,
  canonical_translation_id    INT(10) UNSIGNED,
  stable_id                   VARCHAR(128) DEFAULT NULL,
  version                     SMALLINT UNSIGNED NOT NULL DEFAULT 1,
  created_date                DATETIME NOT NULL DEFAULT '0000-00-00 00:00:00',
  modified_date               DATETIME NOT NULL DEFAULT '0000-00-00 00:00:00',

  PRIMARY KEY (transcript_id),
  KEY seq_region_idx (seq_region_id, seq_region_start),
  KEY gene_index (gene_id),
  KEY xref_id_index (display_xref_id),
  KEY analysis_idx (analysis_id),
  UNIQUE INDEX canonical_translation_idx (canonical_translation_id),
  KEY stable_id_idx (stable_id, version)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;
