CREATE TABLE exon (

  exon_id                     INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  seq_region_id               INT(10) UNSIGNED NOT NULL,
  seq_region_start            INT(10) UNSIGNED NOT NULL,
  seq_region_end              INT(10) UNSIGNED NOT NULL,
  seq_region_strand           TINYINT(2) NOT NULL,

  phase                       TINYINT(2) NOT NULL,
  end_phase                   TINYINT(2) NOT NULL,

  is_current                  BOOLEAN NOT NULL DEFAULT 1,
  is_constitutive             BOOLEAN NOT NULL DEFAULT 0,

  stable_id                   VARCHAR(128) DEFAULT NULL,
  version                     SMALLINT UNSIGNED NOT NULL DEFAULT 1,
  created_date                DATETIME NOT NULL DEFAULT '0000-00-00 00:00:00',
  modified_date               DATETIME NOT NULL DEFAULT '0000-00-00 00:00:00',

  PRIMARY KEY (exon_id),
  KEY seq_region_idx (seq_region_id, seq_region_start),
  KEY stable_id_idx (stable_id, version)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;
