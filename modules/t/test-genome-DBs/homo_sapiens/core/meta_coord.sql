CREATE TABLE `meta_coord` (
  `table_name` varchar(40) collate latin1_bin NOT NULL default '',
  `coord_system_id` int(11) NOT NULL default '0',
  `max_length` int(11) default NULL,
  UNIQUE KEY `table_name` (`table_name`,`coord_system_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

