CREATE TABLE `exon` (
  `exon_id` int(10) unsigned NOT NULL auto_increment,
  `seq_region_id` int(10) unsigned NOT NULL default '0',
  `seq_region_start` int(10) unsigned NOT NULL default '0',
  `seq_region_end` int(10) unsigned NOT NULL default '0',
  `seq_region_strand` tinyint(2) NOT NULL default '0',
  `phase` tinyint(2) NOT NULL default '0',
  `end_phase` tinyint(2) NOT NULL default '0',
  `is_current` BOOLEAN DEFAULT 1,
  PRIMARY KEY  (`exon_id`),
  KEY `seq_region_idx` (`seq_region_id`,`seq_region_start`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

