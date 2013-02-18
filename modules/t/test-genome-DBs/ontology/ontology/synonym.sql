SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `synonym` (
  `synonym_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `term_id` int(10) unsigned NOT NULL,
  `name` text NOT NULL,
  PRIMARY KEY (`synonym_id`),
  UNIQUE KEY `term_synonym_idx` (`term_id`,`synonym_id`),
  KEY `name_idx` (`name`(50))
) ENGINE=MyISAM AUTO_INCREMENT=104283 DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;
