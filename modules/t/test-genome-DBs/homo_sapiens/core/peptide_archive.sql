CREATE TABLE `peptide_archive` (
  `peptide_archive_id` int(11) NOT NULL auto_increment,
  `md5_checksum` varchar(32) default NULL,
  `peptide_seq` mediumtext NOT NULL,
  PRIMARY KEY  (`peptide_archive_id`),
  KEY `checksum` (`md5_checksum`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

