SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `subset` (
  `subset_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(64) NOT NULL,
  `definition` varchar(128) NOT NULL,
  PRIMARY KEY (`subset_id`),
  UNIQUE KEY `name_idx` (`name`)
) ENGINE=MyISAM AUTO_INCREMENT=18 DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;
