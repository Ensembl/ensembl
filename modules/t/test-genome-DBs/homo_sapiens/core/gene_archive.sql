CREATE TABLE `gene_archive` (
  `gene_stable_id` varchar(128) collate latin1_bin NOT NULL default '',
  `gene_version` smallint(6) NOT NULL default '0',
  `transcript_stable_id` varchar(128) collate latin1_bin NOT NULL default '',
  `transcript_version` smallint(6) NOT NULL default '0',
  `translation_stable_id` varchar(128) collate latin1_bin NOT NULL default '',
  `translation_version` smallint(6) NOT NULL default '0',
  `peptide_archive_id` int(11) NOT NULL default '0',
  `mapping_session_id` int(11) NOT NULL default '0',
  KEY `gene_idx` (`gene_stable_id`,`gene_version`),
  KEY `transcript_idx` (`transcript_stable_id`,`transcript_version`),
  KEY `translation_idx` (`translation_stable_id`,`translation_version`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

