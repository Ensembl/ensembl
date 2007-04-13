CREATE TABLE `meta` (
  `meta_id` int(11) NOT NULL auto_increment,
  `meta_key` varchar(40) collate latin1_bin NOT NULL default '',
  `meta_value` varchar(255) collate latin1_bin NOT NULL default '',
  PRIMARY KEY  (`meta_id`),
  KEY `meta_key_index` (`meta_key`),
  KEY `meta_value_index` (`meta_value`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

