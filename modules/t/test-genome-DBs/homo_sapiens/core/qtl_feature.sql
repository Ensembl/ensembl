CREATE TABLE `qtl_feature` (
  `seq_region_id` int(11) NOT NULL default '0',
  `seq_region_start` int(11) NOT NULL default '0',
  `seq_region_end` int(11) NOT NULL default '0',
  `qtl_id` int(11) NOT NULL default '0',
  `analysis_id` int(11) NOT NULL default '0',
  KEY `qtl_id` (`qtl_id`),
  KEY `loc_idx` (`seq_region_id`,`seq_region_start`),
  KEY `analysis_idx` (`analysis_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

