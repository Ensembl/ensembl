SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `relation_type` (
  `relation_type_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(64) NOT NULL,
  PRIMARY KEY (`relation_type_id`),
  UNIQUE KEY `name_idx` (`name`)
) ENGINE=MyISAM AUTO_INCREMENT=85 DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;
