CREATE TABLE `regulatory_feature` (
  `regulatory_feature_id` int(11) NOT NULL auto_increment,
  `name` varchar(255) NOT NULL default '',
  `seq_region_id` int(11) NOT NULL default '0',
  `seq_region_start` int(11) NOT NULL default '0',
  `seq_region_end` int(11) NOT NULL default '0',
  `seq_region_strand` tinyint(4) NOT NULL default '0',
  `analysis_id` int(11) NOT NULL default '0',
  `regulatory_factor_id` int(11) default NULL,
  PRIMARY KEY  (`regulatory_feature_id`),
  KEY `seq_region_idx` (`seq_region_id`,`analysis_id`,`seq_region_start`),
  KEY `seq_region_idx_2` (`seq_region_id`,`seq_region_start`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

