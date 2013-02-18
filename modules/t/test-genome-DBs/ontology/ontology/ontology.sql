SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `ontology` (
  `ontology_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(64) NOT NULL,
  `namespace` varchar(64) NOT NULL,
  PRIMARY KEY (`ontology_id`),
  UNIQUE KEY `name_namespace_idx` (`name`,`namespace`)
) ENGINE=MyISAM AUTO_INCREMENT=8 DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;
