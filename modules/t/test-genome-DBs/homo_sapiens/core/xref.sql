CREATE TABLE `xref` (
  `xref_id` int(10) unsigned NOT NULL auto_increment,
  `external_db_id` int(11) NOT NULL,
  `dbprimary_acc` varchar(40) collate latin1_bin NOT NULL,
  `display_label` varchar(128) collate latin1_bin NOT NULL,
  `version` varchar(10) collate latin1_bin NOT NULL default 0,
  `description` text collate latin1_bin,
  `info_type` enum('NONE','PROJECTION','MISC','DEPENDENT','DIRECT','SEQUENCE_MATCH','INFERRED_PAIR','PROBE','UNMAPPED','CHECKSUM') collate latin1_bin NOT NULL default 'NONE',
  `info_text` varchar(255) collate latin1_bin NOT NULL default '',
  PRIMARY KEY  (`xref_id`),
  UNIQUE KEY `id_index` (`dbprimary_acc`,`external_db_id`,`info_type`,`info_text`),
  KEY `display_index` (`display_label`),
  KEY info_type_idx (info_type)
) ENGINE=MyISAM AUTO_INCREMENT=1000006 DEFAULT CHARSET=latin1 COLLATE=latin1_bin;
