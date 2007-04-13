CREATE TABLE `attrib_type` (
  `attrib_type_id` smallint(5) unsigned NOT NULL auto_increment,
  `code` varchar(15) collate latin1_bin NOT NULL default '',
  `name` varchar(255) collate latin1_bin NOT NULL default '',
  `description` text collate latin1_bin,
  PRIMARY KEY  (`attrib_type_id`),
  UNIQUE KEY `c` (`code`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

