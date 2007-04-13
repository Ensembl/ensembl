CREATE TABLE `misc_attrib` (
  `misc_feature_id` int(10) unsigned NOT NULL default '0',
  `attrib_type_id` smallint(5) unsigned NOT NULL default '0',
  `value` varchar(255) collate latin1_bin NOT NULL default '',
  KEY `type_val_idx` (`attrib_type_id`,`value`),
  KEY `misc_feature_idx` (`misc_feature_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

