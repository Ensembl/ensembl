CREATE TABLE `qtl` (
  `qtl_id` int(10) unsigned NOT NULL auto_increment,
  `trait` varchar(255) collate latin1_bin NOT NULL default '',
  `lod_score` float default NULL,
  `flank_marker_id_1` int(11) default NULL,
  `flank_marker_id_2` int(11) default NULL,
  `peak_marker_id` int(11) default NULL,
  PRIMARY KEY  (`qtl_id`),
  KEY `trait_idx` (`trait`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

