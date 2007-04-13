CREATE TABLE `translation` (
  `translation_id` int(10) unsigned NOT NULL auto_increment,
  `transcript_id` int(10) unsigned NOT NULL default '0',
  `seq_start` int(10) NOT NULL default '0',
  `start_exon_id` int(10) unsigned NOT NULL default '0',
  `seq_end` int(10) NOT NULL default '0',
  `end_exon_id` int(10) unsigned NOT NULL default '0',
  PRIMARY KEY  (`translation_id`),
  KEY `transcript_id` (`transcript_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

