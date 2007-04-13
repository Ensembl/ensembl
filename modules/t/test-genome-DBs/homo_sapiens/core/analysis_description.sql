CREATE TABLE `analysis_description` (
  `analysis_id` int(10) unsigned NOT NULL default '0',
  `description` text collate latin1_bin,
  `display_label` varchar(255) collate latin1_bin default NULL,
  `displayable` BOOLEAN NOT NULL DEFAULT 1,
  `web_data` TEXT,

  KEY `analysis_idx` (`analysis_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

