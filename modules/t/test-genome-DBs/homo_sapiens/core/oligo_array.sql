CREATE TABLE `oligo_array` (
  `oligo_array_id` int(11) NOT NULL auto_increment,
  `parent_array_id` int(11) default NULL,
  `probe_setsize` tinyint(4) NOT NULL default '0',
  `name` varchar(40) collate latin1_bin NOT NULL default '',
  `type` ENUM('AFFY', 'OLIGO'),
  PRIMARY KEY  (`oligo_array_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

