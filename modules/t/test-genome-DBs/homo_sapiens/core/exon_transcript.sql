CREATE TABLE `exon_transcript` (
  `exon_id` int(10) unsigned NOT NULL default '0',
  `transcript_id` int(10) unsigned NOT NULL default '0',
  `rank` int(10) NOT NULL default '0',
  PRIMARY KEY  (`exon_id`,`transcript_id`,`rank`),
  KEY `transcript` (`transcript_id`),
  KEY `exon` (`exon_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

