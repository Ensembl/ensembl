CREATE TABLE `translation_attrib` (
  `translation_id` int(10) unsigned NOT NULL default '0',
  `attrib_type_id` smallint(5) unsigned NOT NULL default '0',
  `value` varchar(255) collate latin1_bin NOT NULL default '',
  KEY `type_val_idx` (`attrib_type_id`,`value`),
  KEY `translation_idx` (`translation_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

