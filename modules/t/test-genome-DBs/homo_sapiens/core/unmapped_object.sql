CREATE TABLE `unmapped_object` (

  `unmapped_object_id`    INT UNSIGNED NOT NULL AUTO_INCREMENT,
  `type`                  ENUM("xref", "cDNA", "Marker") NOT NULL,
  `analysis_id`           INT(10) UNSIGNED NOT NULL,
  `external_db_id`        INT,
  `identifier`            VARCHAR(255) NOT NULL,
  `unmapped_reason_id`    SMALLINT(5) UNSIGNED NOT NULL,
  `query_score`           DOUBLE,
  `target_score`          DOUBLE, 
  `ensembl_id`            INT(10) unsigned default '0',
  `ensembl_object_type`   ENUM('RawContig','Transcript','Gene','Translation') collate latin1_bin default 'RawContig',
  PRIMARY KEY            (`unmapped_object_id`),
  KEY                    id_idx(`identifier`),
  KEY                    anal_idx(`analysis_id`),
  KEY                    anal_exdb_idx(`analysis_id`,`external_db_id`)
   
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;
