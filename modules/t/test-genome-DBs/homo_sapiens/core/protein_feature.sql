CREATE TABLE `protein_feature` (
  `protein_feature_id` int(10) unsigned NOT NULL auto_increment,
  `translation_id` int(11) NOT NULL default '0',
  `seq_start` int(10) NOT NULL default '0',
  `seq_end` int(10) NOT NULL default '0',
  `hit_start` int(10) NOT NULL default '0',
  `hit_end` int(10) NOT NULL default '0',
  `hit_name` varchar(40) collate latin1_bin NOT NULL default '',
  `analysis_id` int(10) unsigned NOT NULL default '0',
  `score` double NOT NULL default '0',
  `evalue` double default NULL,
  `perc_ident` float default NULL,
  PRIMARY KEY  (`protein_feature_id`),
  KEY `translation_id` (`translation_id`),
  KEY `hitname_index` (`hit_name`),
  KEY `analysis_idx` (`analysis_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

