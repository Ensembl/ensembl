CREATE TABLE `protein_align_feature` (
  `protein_align_feature_id` int(10) unsigned NOT NULL auto_increment,
  `seq_region_id` int(10) unsigned NOT NULL default '0',
  `seq_region_start` int(10) unsigned NOT NULL default '0',
  `seq_region_end` int(10) unsigned NOT NULL default '0',
  `seq_region_strand` tinyint(1) NOT NULL default '1',
  `hit_start` int(10) NOT NULL default '0',
  `hit_end` int(10) NOT NULL default '0',
  `hit_name` varchar(40) collate latin1_bin NOT NULL default '',
  `analysis_id` int(10) unsigned NOT NULL default '0',
  `score` double default NULL,
  `evalue` double default NULL,
  `perc_ident` float default NULL,
  `cigar_line` text collate latin1_bin,  
  `external_db_id` smallint(5) unsigned default NULL,
  `hcoverage` double default NULL,
   PRIMARY KEY  (`protein_align_feature_id`),
  KEY `seq_region_idx` (`seq_region_id`,`analysis_id`,`seq_region_start`,`score`),
  KEY `seq_region_idx_2` (`seq_region_id`,`seq_region_start`),
  KEY `hit_idx` (`hit_name`),
  KEY `analysis_idx` (`analysis_id`),
  KEY `external_db_idx` (`external_db_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin MAX_ROWS=100000000 AVG_ROW_LENGTH=80;

