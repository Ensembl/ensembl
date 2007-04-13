CREATE TABLE `interpro` (
  `interpro_ac` varchar(40) collate latin1_bin NOT NULL default '',
  `id` varchar(40) collate latin1_bin NOT NULL default '',
  UNIQUE KEY `interpro_ac` (`interpro_ac`,`id`),
  KEY `id` (`id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

