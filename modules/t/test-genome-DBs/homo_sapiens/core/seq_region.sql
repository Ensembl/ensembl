CREATE TABLE `seq_region` (
  `seq_region_id` int(10) unsigned NOT NULL auto_increment,
  `name` varchar(40) NOT NULL default '',
  `coord_system_id` int(10) NOT NULL default '0',
  `length` int(10) NOT NULL default '0',
  PRIMARY KEY  (`seq_region_id`),
  UNIQUE KEY `coord_system_id` (`coord_system_id`,`name`),
  KEY `name_idx` (`name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

