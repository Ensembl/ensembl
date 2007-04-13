CREATE TABLE `assembly_exception` (
  `assembly_exception_id` int(10) unsigned NOT NULL auto_increment,
  `seq_region_id` int(11) NOT NULL default '0',
  `seq_region_start` int(11) NOT NULL default '0',
  `seq_region_end` int(11) NOT NULL default '0',
  `exc_type` enum('HAP','PAR') collate latin1_bin NOT NULL default 'HAP',
  `exc_seq_region_id` int(11) NOT NULL default '0',
  `exc_seq_region_start` int(11) NOT NULL default '0',
  `exc_seq_region_end` int(11) NOT NULL default '0',
  `ori` int(11) NOT NULL default '0',
  PRIMARY KEY  (`assembly_exception_id`),
  KEY `sr_idx` (`seq_region_id`,`seq_region_start`),
  KEY `ex_idx` (`exc_seq_region_id`,`exc_seq_region_start`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

