CREATE TABLE `identity_xref` (
  `object_xref_id` int(10) unsigned NOT NULL default '0',
  `query_identity` int(5) default NULL,
  `target_identity` int(5) default NULL,
  `hit_start` int(11) default NULL,
  `hit_end` int(11) default NULL,
  `translation_start` int(11) default NULL,
  `translation_end` int(11) default NULL,
  `cigar_line` text collate latin1_bin,
  `score` double default NULL,
  `evalue` double default NULL,
  `analysis_id` int(11) default NULL,
  PRIMARY KEY  (`object_xref_id`),
  KEY `analysis_idx` (`analysis_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

