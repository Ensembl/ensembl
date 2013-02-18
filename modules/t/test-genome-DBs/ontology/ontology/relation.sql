SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `relation` (
  `relation_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `child_term_id` int(10) unsigned NOT NULL,
  `parent_term_id` int(10) unsigned NOT NULL,
  `relation_type_id` int(10) unsigned NOT NULL,
  `intersection_of` tinyint(3) unsigned NOT NULL DEFAULT '0',
  `ontology_id` int(10) unsigned NOT NULL,
  PRIMARY KEY (`relation_id`),
  UNIQUE KEY `child_parent_idx` (`child_term_id`,`parent_term_id`,`relation_type_id`,`intersection_of`,`ontology_id`),
  KEY `parent_idx` (`parent_term_id`)
) ENGINE=MyISAM AUTO_INCREMENT=83305 DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;
