CREATE TABLE `coord_system` (
  `coord_system_id` int(10) unsigned NOT NULL auto_increment,
  `species_id` int(10) unsigned NOT NULL default '1',
  `name` varchar(40) NOT NULL,
  `version` varchar(255) default NULL,
  `rank` int(11) NOT NULL,
  `attrib` set('default_version','sequence_level') default NULL,
  PRIMARY KEY  (`coord_system_id`),
  UNIQUE KEY `rank_idx` (`rank`,`species_id`),
  UNIQUE KEY `name_idx` (`name`,`version`,`species_id`),
  KEY `species_idx` (`species_id`)
) ENGINE=MyISAM AUTO_INCREMENT=8 DEFAULT CHARSET=latin1;
