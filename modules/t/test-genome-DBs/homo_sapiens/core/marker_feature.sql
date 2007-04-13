CREATE TABLE `marker_feature` (
  `marker_feature_id` int(10) unsigned NOT NULL auto_increment,
  `marker_id` int(10) unsigned NOT NULL default '0',
  `seq_region_id` int(10) unsigned NOT NULL default '0',
  `seq_region_start` int(10) unsigned NOT NULL default '0',
  `seq_region_end` int(10) unsigned NOT NULL default '0',
  `analysis_id` int(10) unsigned NOT NULL default '0',
  `map_weight` int(10) unsigned default NULL,
  PRIMARY KEY  (`marker_feature_id`),
  KEY `seq_region_idx` (`seq_region_id`,`seq_region_start`),
  KEY `analysis_idx` (`analysis_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

