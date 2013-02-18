SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `meta` (
  `meta_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `meta_key` varchar(64) NOT NULL,
  `meta_value` varchar(128) DEFAULT NULL,
  PRIMARY KEY (`meta_id`),
  UNIQUE KEY `key_value_idx` (`meta_key`,`meta_value`)
) ENGINE=MyISAM AUTO_INCREMENT=8 DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;
