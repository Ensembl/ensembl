SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `gene_expression` (
  `gene_expression_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `gene_id` int(10) unsigned NOT NULL,
  `tissue_id` int(10) unsigned NOT NULL,
  `value` text NOT NULL,
  `analysis_id` smallint(5) unsigned NOT NULL,
  PRIMARY KEY (`gene_expression_id`),
  UNIQUE KEY `gene_expression_idx` (`gene_id`,`tissue_id`)
) ENGINE=MyISAM AUTO_INCREMENT=52811 DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;
