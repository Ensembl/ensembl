CREATE TABLE `ontology_xref` (
  `object_xref_id` int(10) unsigned NOT NULL default '0',
  `linkage_type` varchar(3) DEFAULT NULL,
  `source_xref_id` int(10) unsigned default NULL,
  UNIQUE KEY `object_xref_id_2` (`object_xref_id`,`source_xref_id`,`linkage_type`),
  KEY `object_xref_id` (`object_xref_id`),
  KEY `source_xref_id` (`source_xref_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
