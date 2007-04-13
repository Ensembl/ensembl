CREATE TABLE `density_feature` (
  `density_feature_id` int(11) NOT NULL auto_increment,
  `density_type_id` int(11) NOT NULL default '0',
  `seq_region_id` int(11) NOT NULL default '0',
  `seq_region_start` int(11) NOT NULL default '0',
  `seq_region_end` int(11) NOT NULL default '0',
  `density_value` float NOT NULL default '0',
  PRIMARY KEY  (`density_feature_id`),
  KEY `seq_region_idx` (`density_type_id`,`seq_region_id`,`seq_region_start`),
  KEY `seq_region_id_idx` (`seq_region_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

