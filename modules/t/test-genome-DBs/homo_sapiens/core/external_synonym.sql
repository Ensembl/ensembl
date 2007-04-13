CREATE TABLE `external_synonym` (
  `xref_id` int(10) unsigned NOT NULL default '0',
  `synonym` varchar(40) collate latin1_bin NOT NULL default '',
  PRIMARY KEY  (`xref_id`,`synonym`),
  KEY `name_index` (`synonym`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

