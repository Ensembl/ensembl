CREATE TABLE `density_type` (
  `density_type_id` int(11) NOT NULL auto_increment,
  `analysis_id` int(11) NOT NULL default '0',
  `block_size` int(11) NOT NULL default '0',
  `region_features` int(11) NOT NULL default '0',
  `value_type` enum('sum','ratio') collate latin1_bin NOT NULL default 'sum',
  PRIMARY KEY  (`density_type_id`),
  UNIQUE KEY `analysis_id` (`analysis_id`,`block_size`,`region_features`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

