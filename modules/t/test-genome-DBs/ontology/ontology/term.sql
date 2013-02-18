SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `term` (
  `term_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `ontology_id` int(10) unsigned NOT NULL,
  `subsets` text,
  `accession` varchar(64) NOT NULL,
  `name` varchar(255) NOT NULL,
  `definition` text,
  `is_root` int(11) DEFAULT NULL,
  PRIMARY KEY (`term_id`),
  UNIQUE KEY `accession_idx` (`accession`),
  UNIQUE KEY `ontology_acc_idx` (`ontology_id`,`accession`),
  KEY `name_idx` (`name`)
) ENGINE=MyISAM AUTO_INCREMENT=44841 DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;
