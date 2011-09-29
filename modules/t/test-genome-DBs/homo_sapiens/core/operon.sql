CREATE TABLE operon (
  operon_id                 INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  seq_region_id             INT(10) UNSIGNED NOT NULL,
  seq_region_start          INT(10) UNSIGNED NOT NULL,
  seq_region_end            INT(10) UNSIGNED NOT NULL,
  seq_region_strand         TINYINT(2) NOT NULL,
  display_label             VARCHAR(255) DEFAULT NULL,
  analysis_id               SMALLINT UNSIGNED NOT NULL,
  stable_id                 VARCHAR(128) DEFAULT NULL,
  version                   SMALLINT UNSIGNED NOT NULL DEFAULT 1,
  created_date              DATETIME NOT NULL DEFAULT '0000-00-00 00:00:00',
  modified_date             DATETIME NOT NULL DEFAULT '0000-00-00 00:00:00',

  PRIMARY KEY (operon_id),
  KEY seq_region_idx (seq_region_id, seq_region_start),
  KEY name_idx (display_label),
  KEY stable_id_idx (stable_id, version)
) COLLATE=latin1_swedish_ci ENGINE=MyISAM;
