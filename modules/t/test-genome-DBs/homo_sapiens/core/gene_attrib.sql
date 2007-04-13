CREATE TABLE `gene_attrib` (
  `gene_id` int(10) unsigned NOT NULL default '0',
  `attrib_type_id` smallint(5) unsigned NOT NULL default '0',
  `value` varchar(255) NOT NULL default '',
  KEY `type_val_idx` (`attrib_type_id`,`value`),
  KEY `gene_idx` (`gene_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;
