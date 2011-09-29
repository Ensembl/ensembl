CREATE TABLE translation (

  translation_id              INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  transcript_id               INT(10) UNSIGNED NOT NULL,
  seq_start                   INT(10) NOT NULL,       # relative to exon start
  start_exon_id               INT(10) UNSIGNED NOT NULL,
  seq_end                     INT(10) NOT NULL,       # relative to exon start
  end_exon_id                 INT(10) UNSIGNED NOT NULL,
  stable_id                   VARCHAR(128) DEFAULT NULL,
  version                     SMALLINT UNSIGNED NOT NULL DEFAULT 1,
  created_date                DATETIME NOT NULL DEFAULT '0000-00-00 00:00:00',
  modified_date               DATETIME NOT NULL DEFAULT '0000-00-00 00:00:00',

  PRIMARY KEY (translation_id),
  KEY transcript_idx (transcript_id),
  KEY stable_id_idx (stable_id, version)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;
