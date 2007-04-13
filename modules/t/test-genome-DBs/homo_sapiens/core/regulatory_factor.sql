CREATE TABLE `regulatory_factor` (
  `regulatory_factor_id` int(11) NOT NULL auto_increment,
  `name` varchar(255) NOT NULL default '',
  `type` enum('miRNA_target','transcription_factor','transcription_factor_complex') default NULL,
  PRIMARY KEY  (`regulatory_factor_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

