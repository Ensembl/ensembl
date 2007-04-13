CREATE TABLE `repeat_consensus` (
  `repeat_consensus_id` int(10) unsigned NOT NULL auto_increment,
  `repeat_name` varchar(255) collate latin1_bin NOT NULL default '',
  `repeat_class` varchar(100) collate latin1_bin NOT NULL default '',
  `repeat_type` varchar(40) collate latin1_bin NOT NULL default '',
  `repeat_consensus` text collate latin1_bin,
  PRIMARY KEY  (`repeat_consensus_id`),
  KEY `name` (`repeat_name`),
  KEY `class` (`repeat_class`),
  KEY `consensus` (`repeat_consensus`(10)),
  KEY `type` (`repeat_type`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

