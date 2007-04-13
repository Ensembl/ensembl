CREATE TABLE `regulatory_search_region` (
  `regulatory_search_region_id` int(11) NOT NULL auto_increment,
  `name` varchar(255) NOT NULL default '',
  `seq_region_id` int(11) NOT NULL default '0',
  `seq_region_start` int(11) NOT NULL default '0',
  `seq_region_end` int(11) NOT NULL default '0',
  `seq_region_strand` tinyint(4) NOT NULL default '0',
  `ensembl_object_type` enum('Transcript','Translation','Gene') NOT NULL default 'Transcript',
  `ensembl_object_id` int(11) default NULL,
  `analysis_id` int(11) NOT NULL default '0',
  PRIMARY KEY  (`regulatory_search_region_id`),
  KEY `rsr_idx` (`regulatory_search_region_id`),
  KEY `ensembl_object_idx` (`ensembl_object_type`,`ensembl_object_id`),
  KEY `seq_region_idx` (`seq_region_id`,`seq_region_start`),
  KEY `seq_region_idx_2` (`seq_region_id`,`seq_region_start`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

