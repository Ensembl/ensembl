CREATE TABLE `alt_id` (
  `alt_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `term_id` int(10) unsigned NOT NULL,
  `accession` varchar(64) NOT NULL,
  PRIMARY KEY (`alt_id`),
  UNIQUE KEY `term_alt_idx` (`term_id`,`alt_id`),
  KEY `ix_alt_id_accession` (`accession`)
) ENGINE=MyISAM AUTO_INCREMENT=7 DEFAULT CHARSET=latin1;

CREATE TABLE `aux_GO_Cross_product_review_map` (
  `term_id` int(10) unsigned NOT NULL,
  `subset_term_id` int(10) unsigned NOT NULL,
  `distance` tinyint(3) unsigned NOT NULL,
  UNIQUE KEY `map_idx` (`term_id`,`subset_term_id`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8;

CREATE TABLE `aux_GO_goslim_aspergillus_map` (
  `term_id` int(10) unsigned NOT NULL,
  `subset_term_id` int(10) unsigned NOT NULL,
  `distance` tinyint(3) unsigned NOT NULL,
  UNIQUE KEY `map_idx` (`term_id`,`subset_term_id`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8;

CREATE TABLE `aux_GO_goslim_candida_map` (
  `term_id` int(10) unsigned NOT NULL,
  `subset_term_id` int(10) unsigned NOT NULL,
  `distance` tinyint(3) unsigned NOT NULL,
  UNIQUE KEY `map_idx` (`term_id`,`subset_term_id`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8;

CREATE TABLE `aux_GO_goslim_generic_map` (
  `term_id` int(10) unsigned NOT NULL,
  `subset_term_id` int(10) unsigned NOT NULL,
  `distance` tinyint(3) unsigned NOT NULL,
  UNIQUE KEY `map_idx` (`term_id`,`subset_term_id`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8;

CREATE TABLE `aux_GO_goslim_metagenomics_map` (
  `term_id` int(10) unsigned NOT NULL,
  `subset_term_id` int(10) unsigned NOT NULL,
  `distance` tinyint(3) unsigned NOT NULL,
  UNIQUE KEY `map_idx` (`term_id`,`subset_term_id`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8;

CREATE TABLE `aux_GO_goslim_pir_map` (
  `term_id` int(10) unsigned NOT NULL,
  `subset_term_id` int(10) unsigned NOT NULL,
  `distance` tinyint(3) unsigned NOT NULL,
  UNIQUE KEY `map_idx` (`term_id`,`subset_term_id`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8;

CREATE TABLE `aux_GO_goslim_plant_map` (
  `term_id` int(10) unsigned NOT NULL,
  `subset_term_id` int(10) unsigned NOT NULL,
  `distance` tinyint(3) unsigned NOT NULL,
  UNIQUE KEY `map_idx` (`term_id`,`subset_term_id`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8;

CREATE TABLE `aux_GO_goslim_pombe_map` (
  `term_id` int(10) unsigned NOT NULL,
  `subset_term_id` int(10) unsigned NOT NULL,
  `distance` tinyint(3) unsigned NOT NULL,
  UNIQUE KEY `map_idx` (`term_id`,`subset_term_id`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8;

CREATE TABLE `aux_GO_goslim_yeast_map` (
  `term_id` int(10) unsigned NOT NULL,
  `subset_term_id` int(10) unsigned NOT NULL,
  `distance` tinyint(3) unsigned NOT NULL,
  UNIQUE KEY `map_idx` (`term_id`,`subset_term_id`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8;

CREATE TABLE `aux_GO_gosubset_prok_map` (
  `term_id` int(10) unsigned NOT NULL,
  `subset_term_id` int(10) unsigned NOT NULL,
  `distance` tinyint(3) unsigned NOT NULL,
  UNIQUE KEY `map_idx` (`term_id`,`subset_term_id`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8;

CREATE TABLE `aux_GO_high_level_annotation_qc_map` (
  `term_id` int(10) unsigned NOT NULL,
  `subset_term_id` int(10) unsigned NOT NULL,
  `distance` tinyint(3) unsigned NOT NULL,
  UNIQUE KEY `map_idx` (`term_id`,`subset_term_id`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8;

CREATE TABLE `aux_GO_mf_needs_review_map` (
  `term_id` int(10) unsigned NOT NULL,
  `subset_term_id` int(10) unsigned NOT NULL,
  `distance` tinyint(3) unsigned NOT NULL,
  UNIQUE KEY `map_idx` (`term_id`,`subset_term_id`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8;

CREATE TABLE `aux_GO_virus_checked_map` (
  `term_id` int(10) unsigned NOT NULL,
  `subset_term_id` int(10) unsigned NOT NULL,
  `distance` tinyint(3) unsigned NOT NULL,
  UNIQUE KEY `map_idx` (`term_id`,`subset_term_id`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8;

CREATE TABLE `aux_SO_DBVAR_map` (
  `term_id` int(10) unsigned NOT NULL,
  `subset_term_id` int(10) unsigned NOT NULL,
  `distance` tinyint(3) unsigned NOT NULL,
  UNIQUE KEY `map_idx` (`term_id`,`subset_term_id`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8;

CREATE TABLE `aux_SO_SOFA_map` (
  `term_id` int(10) unsigned NOT NULL,
  `subset_term_id` int(10) unsigned NOT NULL,
  `distance` tinyint(3) unsigned NOT NULL,
  UNIQUE KEY `map_idx` (`term_id`,`subset_term_id`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8;

CREATE TABLE `aux_SO_biosapiens_map` (
  `term_id` int(10) unsigned NOT NULL,
  `subset_term_id` int(10) unsigned NOT NULL,
  `distance` tinyint(3) unsigned NOT NULL,
  UNIQUE KEY `map_idx` (`term_id`,`subset_term_id`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8;

CREATE TABLE `closure` (
  `closure_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `child_term_id` int(10) unsigned NOT NULL,
  `parent_term_id` int(10) unsigned NOT NULL,
  `subparent_term_id` int(10) unsigned DEFAULT NULL,
  `distance` tinyint(3) unsigned NOT NULL,
  `ontology_id` int(10) unsigned NOT NULL,
  `confident_relationship` tinyint(1) NOT NULL,
  PRIMARY KEY (`closure_id`),
  UNIQUE KEY `closure_child_parent_idx` (`child_term_id`,`parent_term_id`,`subparent_term_id`,`ontology_id`),
  KEY `ix_closure_subparent_term_id` (`subparent_term_id`),
  KEY `ix_closure_ontology_id` (`ontology_id`),
  KEY `parent_subparent_idx` (`parent_term_id`,`subparent_term_id`)
) ENGINE=MyISAM AUTO_INCREMENT=1453438 DEFAULT CHARSET=latin1;

CREATE TABLE `meta` (
  `meta_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `meta_key` varchar(64) COLLATE utf8_unicode_ci NOT NULL,
  `meta_value` varchar(128) COLLATE utf8_unicode_ci DEFAULT NULL,
  `species_id` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`meta_id`),
  UNIQUE KEY `key_value_idx` (`meta_key`,`meta_value`)
) ENGINE=MyISAM AUTO_INCREMENT=46 DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;

CREATE TABLE `ontology` (
  `ontology_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(64) COLLATE utf8_unicode_ci NOT NULL,
  `namespace` varchar(64) COLLATE utf8_unicode_ci NOT NULL,
  `data_version` varchar(64) COLLATE utf8_unicode_ci DEFAULT NULL,
  `title` varchar(255) COLLATE utf8_unicode_ci DEFAULT NULL,
  PRIMARY KEY (`ontology_id`),
  UNIQUE KEY `ontology_name_namespace_idx` (`name`,`namespace`)
) ENGINE=MyISAM AUTO_INCREMENT=8 DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;

CREATE TABLE `relation` (
  `relation_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `child_term_id` int(10) unsigned NOT NULL,
  `parent_term_id` int(10) unsigned NOT NULL,
  `relation_type_id` int(10) unsigned NOT NULL,
  `intersection_of` tinyint(1) NOT NULL,
  `ontology_id` int(10) unsigned NOT NULL,
  PRIMARY KEY (`relation_id`),
  UNIQUE KEY `child_parent_idx` (`child_term_id`,`parent_term_id`,`relation_type_id`,`intersection_of`,`ontology_id`),
  KEY `ix_relation_parent_term_id` (`parent_term_id`),
  KEY `ix_relation_relation_type_id` (`relation_type_id`),
  KEY `ix_relation_ontology_id` (`ontology_id`)
) ENGINE=MyISAM AUTO_INCREMENT=68750 DEFAULT CHARSET=latin1;

CREATE TABLE `relation_type` (
  `relation_type_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(64) NOT NULL,
  PRIMARY KEY (`relation_type_id`),
  UNIQUE KEY `name` (`name`)
) ENGINE=MyISAM AUTO_INCREMENT=85 DEFAULT CHARSET=latin1;

CREATE TABLE `subset` (
  `subset_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(64) COLLATE utf8_unicode_ci NOT NULL,
  `definition` varchar(1023) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  PRIMARY KEY (`subset_id`),
  UNIQUE KEY `name` (`name`)
) ENGINE=MyISAM AUTO_INCREMENT=18 DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;

CREATE TABLE `synonym` (
  `synonym_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `term_id` int(10) unsigned NOT NULL,
  `name` text CHARACTER SET utf8 NOT NULL,
  `type` enum('EXACT','BROAD','NARROW','RELATED') COLLATE utf8_unicode_ci DEFAULT NULL,
  `dbxref` varchar(500) COLLATE utf8_unicode_ci DEFAULT NULL,
  PRIMARY KEY (`synonym_id`),
  UNIQUE KEY `synonym_term_idx` (`term_id`,`synonym_id`),
  KEY `synonym_name_idx` (`name`(100))
) ENGINE=MyISAM AUTO_INCREMENT=104283 DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;

CREATE TABLE `term` (
  `term_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `ontology_id` int(10) unsigned NOT NULL,
  `subsets` text COLLATE utf8_unicode_ci DEFAULT NULL,
  `accession` varchar(64) COLLATE utf8_unicode_ci NOT NULL,
  `name` varchar(255) CHARACTER SET utf8 NOT NULL,
  `definition` text CHARACTER SET utf8 DEFAULT NULL,
  `is_root` int(11) NOT NULL DEFAULT 0,
  `is_obsolete` int(11) NOT NULL DEFAULT 0,
  `iri` text COLLATE utf8_unicode_ci DEFAULT NULL,
  PRIMARY KEY (`term_id`),
  UNIQUE KEY `accession` (`accession`),
  UNIQUE KEY `term_ontology_acc_idx` (`ontology_id`,`accession`),
  KEY `term_name_idx` (`name`(100))
) ENGINE=MyISAM AUTO_INCREMENT=45001 DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;

