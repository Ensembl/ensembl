CREATE TABLE `analysis` (
  `analysis_id` int(10) unsigned NOT NULL auto_increment,
  `created` datetime NOT NULL default '0000-00-00 00:00:00',
  `logic_name` varchar(40) collate latin1_bin NOT NULL default '',
  `db` varchar(120) collate latin1_bin default NULL,
  `db_version` varchar(40) collate latin1_bin default NULL,
  `db_file` varchar(120) collate latin1_bin default NULL,
  `program` varchar(80) collate latin1_bin default NULL,
  `program_version` varchar(40) collate latin1_bin default NULL,
  `program_file` varchar(80) collate latin1_bin default NULL,
  `parameters` varchar(255) collate latin1_bin default NULL,
  `module` varchar(80) collate latin1_bin default NULL,
  `module_version` varchar(40) collate latin1_bin default NULL,
  `gff_source` varchar(40) collate latin1_bin default NULL,
  `gff_feature` varchar(40) collate latin1_bin default NULL,
  PRIMARY KEY  (`analysis_id`),
  UNIQUE KEY `logic_name` (`logic_name`),
  KEY `logic_name_idx` (`logic_name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

