CREATE TABLE `misc_set` (
  `misc_set_id` smallint(5) unsigned NOT NULL auto_increment,
  `code` varchar(25) collate latin1_bin NOT NULL default '',
  `name` varchar(255) collate latin1_bin NOT NULL default '',
  `description` text collate latin1_bin NOT NULL,
  `max_length` int(10) unsigned NOT NULL default '0',
  PRIMARY KEY  (`misc_set_id`),
  UNIQUE KEY `c` (`code`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

