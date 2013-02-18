SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `aux_GO_goslim_generic_map` (
  `term_id` int(10) unsigned NOT NULL,
  `subset_term_id` int(10) unsigned NOT NULL,
  `distance` tinyint(3) unsigned NOT NULL,
  UNIQUE KEY `map_idx` (`term_id`,`subset_term_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;
