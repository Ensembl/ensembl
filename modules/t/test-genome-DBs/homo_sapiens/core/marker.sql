CREATE TABLE `marker` (
  `marker_id` int(10) unsigned NOT NULL auto_increment,
  `display_marker_synonym_id` int(10) unsigned default NULL,
  `left_primer` varchar(100) collate latin1_bin NOT NULL default '',
  `right_primer` varchar(100) collate latin1_bin NOT NULL default '',
  `min_primer_dist` int(10) unsigned NOT NULL default '0',
  `max_primer_dist` int(10) unsigned NOT NULL default '0',
  `priority` int(11) default NULL,
  `type` enum('est','microsatellite') collate latin1_bin default NULL,
  PRIMARY KEY  (`marker_id`),
  KEY `marker_idx` (`marker_id`,`priority`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

