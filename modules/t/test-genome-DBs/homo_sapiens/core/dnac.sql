CREATE TABLE `dnac` (
  `seq_region_id` int(10) unsigned NOT NULL default '0',
  `sequence` mediumblob NOT NULL,
  `n_line` text collate latin1_bin,
  PRIMARY KEY  (`seq_region_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin MAX_ROWS=750000 AVG_ROW_LENGTH=19000;

