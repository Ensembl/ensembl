CREATE TABLE `oligo_feature` (
  `oligo_feature_id` int(11) NOT NULL auto_increment,
  `seq_region_id` int(10) unsigned NOT NULL default '0',
  `seq_region_start` int(11) NOT NULL default '0',
  `seq_region_end` int(11) NOT NULL default '0',
  `seq_region_strand` tinyint(4) NOT NULL default '0',
  `mismatches` tinyint(4) default NULL,
  `oligo_probe_id` int(11) NOT NULL default '0',
  `analysis_id` int(11) NOT NULL default '0',
  PRIMARY KEY  (`oligo_feature_id`),
  KEY `seq_region_idx` (`seq_region_id`,`seq_region_start`),
  KEY `probe_idx` (`oligo_probe_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

