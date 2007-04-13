CREATE TABLE `regulatory_feature_object` (
  `regulatory_feature_id` int(11) NOT NULL default '0',
  `ensembl_object_type` enum('Transcript','Translation','Gene') NOT NULL default 'Transcript',
  `ensembl_object_id` int(11) NOT NULL default '0',
  `influence` enum('positive','negative','mixed','unknown') default NULL,
  `evidence` varchar(255) default NULL,
  KEY `regulatory_feature_idx` (`regulatory_feature_id`),
  KEY `ensembl_object_idx` (`ensembl_object_type`,`ensembl_object_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

