CREATE TABLE `gene` (
  `gene_id` int(10) unsigned NOT NULL auto_increment,
  `biotype` varchar(40) collate latin1_bin NOT NULL default '',
  `analysis_id` int(11) default NULL,
  `seq_region_id` int(10) unsigned NOT NULL default '0',
  `seq_region_start` int(10) unsigned NOT NULL default '0',
  `seq_region_end` int(10) unsigned NOT NULL default '0',
  `seq_region_strand` tinyint(2) NOT NULL default '0',
  `display_xref_id` int(10) unsigned default NULL,
  `source` varchar(20) collate latin1_bin NOT NULL default '',
  `status` enum('KNOWN','NOVEL','PUTATIVE','PREDICTED','KNOWN_BY_PROJECTION') collate latin1_bin default NULL,
  `description` text collate latin1_bin,
  `is_current` BOOLEAN DEFAULT 1,
  canonical_transcript_id     INT(10) UNSIGNED NOT NULL,
  canonical_annotation        VARCHAR(255) DEFAULT NULL,

  PRIMARY KEY  (`gene_id`),
  KEY `seq_region_idx` (`seq_region_id`,`seq_region_start`),
  KEY `xref_id_index` (`display_xref_id`),
  KEY `analysis_idx` (`analysis_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

