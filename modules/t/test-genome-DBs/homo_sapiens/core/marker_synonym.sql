CREATE TABLE `marker_synonym` (
  `marker_synonym_id` int(10) unsigned NOT NULL auto_increment,
  `marker_id` int(10) unsigned NOT NULL default '0',
  `source` varchar(20) collate latin1_bin default NULL,
  `name` varchar(30) collate latin1_bin default NULL,
  PRIMARY KEY  (`marker_synonym_id`),
  KEY `marker_synonym_idx` (`marker_synonym_id`,`name`),
  KEY `marker_idx` (`marker_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

