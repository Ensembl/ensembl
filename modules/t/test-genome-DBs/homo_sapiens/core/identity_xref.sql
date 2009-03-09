CREATE TABLE `identity_xref` (
  `object_xref_id` int(10) unsigned NOT NULL default '0',
  `xref_identity` int(5) default NULL,
  `ensembl_identity` int(5) default NULL,
  `xref_start` int(11) default NULL,
  `xref_end` int(11) default NULL,
  `ensembl_start` int(11) default NULL,
  `ensembl_end` int(11) default NULL,
  `cigar_line` text collate latin1_bin,
  `score` double default NULL,
  `evalue` double default NULL,
  PRIMARY KEY  (`object_xref_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

