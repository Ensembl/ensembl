CREATE TABLE `oligo_probe` (
  `oligo_probe_id` int(11) NOT NULL auto_increment,
  `oligo_array_id` int(11) NOT NULL default '0',
  `probeset` varchar(40) collate latin1_bin default NULL,
  `name` varchar(20) collate latin1_bin default NULL,
  `description` text collate latin1_bin default NULL,
  `length` smallint(5) NOT NULL default '0',
  PRIMARY KEY  (`oligo_probe_id`,`oligo_array_id`),
  KEY `probeset_idx` (`probeset`),
  KEY `array_idx` (`oligo_array_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

