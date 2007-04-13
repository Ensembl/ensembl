CREATE TABLE `repeat_feature` (
  `repeat_feature_id` int(10) unsigned NOT NULL auto_increment,
  `seq_region_id` int(10) unsigned NOT NULL default '0',
  `seq_region_start` int(10) unsigned NOT NULL default '0',
  `seq_region_end` int(10) unsigned NOT NULL default '0',
  `seq_region_strand` tinyint(1) NOT NULL default '1',
  `repeat_start` int(10) NOT NULL default '0',
  `repeat_end` int(10) NOT NULL default '0',
  `repeat_consensus_id` int(10) unsigned NOT NULL default '0',
  `analysis_id` int(10) unsigned NOT NULL default '0',
  `score` double default NULL,
  PRIMARY KEY  (`repeat_feature_id`),
  KEY `seq_region_idx` (`seq_region_id`,`seq_region_start`),
  KEY `repeat_idx` (`repeat_consensus_id`),
  KEY `analysis_idx` (`analysis_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin MAX_ROWS=100000000 AVG_ROW_LENGTH=80;

