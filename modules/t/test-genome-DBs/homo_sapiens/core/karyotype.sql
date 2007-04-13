CREATE TABLE `karyotype` (
  `karyotype_id` int(10) unsigned NOT NULL auto_increment,
  `seq_region_id` int(10) unsigned NOT NULL default '0',
  `seq_region_start` int(10) NOT NULL default '0',
  `seq_region_end` int(10) NOT NULL default '0',
  `band` varchar(40) collate latin1_bin NOT NULL default '',
  `stain` varchar(40) collate latin1_bin NOT NULL default '',
  PRIMARY KEY  (`karyotype_id`),
  KEY `region_band_idx` (`seq_region_id`,`band`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

