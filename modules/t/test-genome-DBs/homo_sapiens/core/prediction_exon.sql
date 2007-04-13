CREATE TABLE `prediction_exon` (
  `prediction_exon_id` int(10) unsigned NOT NULL auto_increment,
  `prediction_transcript_id` int(10) unsigned NOT NULL default '0',
  `exon_rank` smallint(5) unsigned NOT NULL default '0',
  `seq_region_id` int(10) unsigned NOT NULL default '0',
  `seq_region_start` int(10) unsigned NOT NULL default '0',
  `seq_region_end` int(10) unsigned NOT NULL default '0',
  `seq_region_strand` tinyint(4) NOT NULL default '0',
  `start_phase` tinyint(4) NOT NULL default '0',
  `score` double default NULL,
  `p_value` double default NULL,
  PRIMARY KEY  (`prediction_exon_id`),
  KEY `prediction_transcript_id` (`prediction_transcript_id`),
  KEY `seq_region_id` (`seq_region_id`,`seq_region_start`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

