CREATE TABLE `simple_feature` (
  `simple_feature_id` int(10) unsigned NOT NULL auto_increment,
  `seq_region_id` int(10) unsigned NOT NULL default '0',
  `seq_region_start` int(10) unsigned NOT NULL default '0',
  `seq_region_end` int(10) unsigned NOT NULL default '0',
  `seq_region_strand` tinyint(1) NOT NULL default '0',
  `display_label` varchar(40) collate latin1_bin NOT NULL default '',
  `analysis_id` int(10) unsigned NOT NULL default '0',
  `score` double default NULL,
  PRIMARY KEY  (`simple_feature_id`),
  KEY `seq_region_idx` (`seq_region_id`,`seq_region_start`),
  KEY `analysis_idx` (`analysis_id`),
  KEY `hit_idx` (`display_label`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin MAX_ROWS=100000000 AVG_ROW_LENGTH=80;

