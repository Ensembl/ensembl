CREATE TABLE `object_xref` (
  `object_xref_id` int(11) NOT NULL auto_increment,
  `ensembl_id` int(10) unsigned NOT NULL default '0',
  `ensembl_object_type` enum('RawContig','Transcript','Gene','Translation','regulatory_factor','regulatory_feature') collate latin1_bin NOT NULL default 'RawContig',
  `xref_id` int(10) unsigned NOT NULL default '0',
  `linkage_annotation` VARCHAR(255) DEFAULT NULL,
  analysis_id                 SMALLINT UNSIGNED NOT NULL,
  UNIQUE KEY `ensembl_object_type` (`ensembl_object_type`,`ensembl_id`,`xref_id`),
  KEY `oxref_idx` (`object_xref_id`,`xref_id`,`ensembl_object_type`,`ensembl_id`),
  KEY `xref_idx` (`xref_id`,`ensembl_object_type`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

