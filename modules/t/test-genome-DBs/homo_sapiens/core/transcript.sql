CREATE TABLE `transcript` (
  `transcript_id` int(10) unsigned NOT NULL auto_increment,
  `gene_id` int(10) unsigned NOT NULL default '0',
  `analysis_id` INT(10) UNSIGNED NOT NULL,
  `seq_region_id` int(10) unsigned NOT NULL default '0',
  `seq_region_start` int(10) unsigned NOT NULL default '0',
  `seq_region_end` int(10) unsigned NOT NULL default '0',
  `seq_region_strand` tinyint(2) NOT NULL default '0',
  `display_xref_id` int(10) unsigned default NULL,
  `biotype` varchar(40) NOT NULL default '',
  `status` enum('KNOWN','NOVEL','PUTATIVE','PREDICTED','KNOWN_BY_PROJECTION') collate latin1_bin default NULL,
  `description` text collate latin1_bin,
  `is_current` BOOLEAN DEFAULT 1,
  PRIMARY KEY  (`transcript_id`),
  KEY `seq_region_idx` (`seq_region_id`,`seq_region_start`),
  KEY `gene_index` (`gene_id`),
  KEY `xref_id_index` (`display_xref_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

