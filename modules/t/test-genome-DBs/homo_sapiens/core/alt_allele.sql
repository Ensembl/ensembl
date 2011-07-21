CREATE TABLE `alt_allele` (
  `alt_allele_id` int(11) NOT NULL auto_increment,
  `gene_id` int(11) NOT NULL default '0',
  `is_ref`  BOOLEAN NOT NULL DEFAULT '0',
  UNIQUE KEY `gene_idx` (`gene_id`),
  UNIQUE KEY `allele_idx` (`alt_allele_id`,`gene_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

