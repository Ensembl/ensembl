CREATE TABLE `regulatory_factor_coding` (
  `regulatory_factor_id` int(11) NOT NULL default '0',
  `transcript_id` int(11) default NULL,
  `gene_id` int(11) default NULL,
  KEY `transcript_idx` (`transcript_id`),
  KEY `gene_idx` (`gene_id`),
  KEY `regulatory_factor_idx` (`regulatory_factor_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

