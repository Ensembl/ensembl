CREATE TABLE `transcript_supporting_feature` (
  `transcript_id` int(11) NOT NULL default '0',
  `feature_type` enum('dna_align_feature','protein_align_feature') collate latin1_bin default NULL,
  `feature_id` int(11) NOT NULL default '0',
  UNIQUE KEY `all_idx` (`transcript_id`,`feature_type`,`feature_id`),
  KEY `feature_idx` (`feature_type`,`feature_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin MAX_ROWS=100000000 AVG_ROW_LENGTH=80;

