CREATE TABLE `qtl_synonym` (
  `qtl_synonym_id` int(10) unsigned NOT NULL auto_increment,
  `qtl_id` int(10) unsigned NOT NULL default '0',
  `source_database` enum('rat genome database','ratmap') collate latin1_bin NOT NULL default 'rat genome database',
  `source_primary_id` varchar(255) collate latin1_bin NOT NULL default '',
  PRIMARY KEY  (`qtl_synonym_id`),
  KEY `qtl_idx` (`qtl_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

