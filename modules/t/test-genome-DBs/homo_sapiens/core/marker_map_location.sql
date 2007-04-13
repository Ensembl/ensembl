CREATE TABLE `marker_map_location` (
  `marker_id` int(10) unsigned NOT NULL default '0',
  `map_id` int(10) unsigned NOT NULL default '0',
  `chromosome_name` varchar(15) collate latin1_bin NOT NULL default '',
  `marker_synonym_id` int(10) unsigned NOT NULL default '0',
  `position` varchar(15) collate latin1_bin NOT NULL default '',
  `lod_score` double default NULL,
  PRIMARY KEY  (`marker_id`,`map_id`),
  KEY `map_idx` (`map_id`,`chromosome_name`,`position`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

