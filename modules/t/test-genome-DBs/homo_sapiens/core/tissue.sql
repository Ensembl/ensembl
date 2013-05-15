SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `tissue` (
  `tissue_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(255) DEFAULT NULL,
  `description` text,
  `is_a` varchar(255) DEFAULT NULL,
  `ontology` varchar(64) NOT NULL,
  PRIMARY KEY (`tissue_id`),
  UNIQUE KEY `name_idx` (`name`)
) ENGINE=MyISAM AUTO_INCREMENT=3 DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;
