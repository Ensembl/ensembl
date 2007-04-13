CREATE TABLE `gene_stable_id` (
  `gene_id` int(10) unsigned NOT NULL default '0',
  `stable_id` varchar(128) collate latin1_bin NOT NULL default '',
  `version` int(10) default NULL,
  `created_date` datetime NOT NULL default '0000-00-00 00:00:00',
  `modified_date` datetime NOT NULL default '0000-00-00 00:00:00',
  PRIMARY KEY  (`gene_id`),
  KEY `stable_id` (`stable_id`,`version`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

