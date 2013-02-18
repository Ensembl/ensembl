SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `closure` (
  `closure_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `child_term_id` int(10) unsigned NOT NULL,
  `parent_term_id` int(10) unsigned NOT NULL,
  `subparent_term_id` int(10) unsigned DEFAULT NULL,
  `distance` tinyint(3) unsigned NOT NULL,
  `ontology_id` int(10) unsigned NOT NULL,
  PRIMARY KEY (`closure_id`),
  UNIQUE KEY `child_parent_idx` (`child_term_id`,`parent_term_id`,`subparent_term_id`,`ontology_id`),
  KEY `parent_subparent_idx` (`parent_term_id`,`subparent_term_id`)
) ENGINE=MyISAM AUTO_INCREMENT=1469494 DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;
