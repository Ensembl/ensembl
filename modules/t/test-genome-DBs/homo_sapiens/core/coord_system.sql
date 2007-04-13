CREATE TABLE `coord_system` (
  `coord_system_id` int(11) NOT NULL auto_increment,
  `name` varchar(40) collate latin1_bin NOT NULL default '',
  `version` varchar(40) collate latin1_bin default NULL,
  `rank` int(11) NOT NULL default '0',
  `attrib` set('default_version','sequence_level') collate latin1_bin default NULL,
  PRIMARY KEY  (`coord_system_id`),
  UNIQUE KEY `rank` (`rank`),
  UNIQUE KEY `name` (`name`,`version`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

