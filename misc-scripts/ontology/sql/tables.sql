-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2018] EMBL-European Bioinformatics Institute
--
-- Licensed under the Apache License, Version 2.0 (the "License");
-- you may not use this file except in compliance with the License.
-- You may obtain a copy of the License at
--
--      http://www.apache.org/licenses/LICENSE-2.0
--
-- Unless required by applicable law or agreed to in writing, software
-- distributed under the License is distributed on an "AS IS" BASIS,
-- WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
-- See the License for the specific language governing permissions and
-- limitations under the License.


-- ---------------------------------------------------------------------
-- The schema for the ensembl_ontology_NN database.
-- ---------------------------------------------------------------------

CREATE TABLE `meta` (
  `meta_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `meta_key` varchar(64) NOT NULL,
  `meta_value` varchar(128) DEFAULT NULL,
  `species_id` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`meta_id`),
  UNIQUE KEY `key_value_idx` (`meta_key`,`meta_value`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;

# Add schema type and schema version to the meta table
INSERT INTO meta (meta_key, meta_value) VALUES
  ('schema_type', 'ontology'),
  ('schema_version', '96');

# Patches included in this schema file
INSERT INTO meta (meta_key, meta_value)
  VALUES ('patch', 'patch_95_96_a.sql|schema_version');


CREATE TABLE `ontology` (
  `ontology_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(64) NOT NULL,
  `namespace` varchar(64) NOT NULL,
  `data_version` varchar(64) DEFAULT NULL,
  `title` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`ontology_id`),
  UNIQUE KEY `ontology_name_namespace_idx` (`name`,`namespace`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;


CREATE TABLE `subset` (
  `subset_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(64) NOT NULL,
  `definition` varchar(511) NOT NULL DEFAULT '',
  PRIMARY KEY (`subset_id`),
  UNIQUE KEY `name` (`name`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;

CREATE TABLE `term` (
  `term_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `ontology_id` int(10) unsigned NOT NULL,
  `subsets` text,
  `accession` varchar(64) NOT NULL,
  `name` varchar(255) COLLATE utf8_general_ci NOT NULL,
  `definition` text COLLATE utf8_general_ci,
  `is_root` int(11) NOT NULL DEFAULT 0,
  `is_obsolete` int(11) NOT NULL DEFAULT 0,
  `iri` text,
  PRIMARY KEY (`term_id`),
  UNIQUE KEY `accession` (`accession`),
  UNIQUE KEY `term_ontology_acc_idx` (`ontology_id`,`accession`),
  KEY `term_name_idx` (`name`(100))
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;

CREATE TABLE `synonym` (
  `synonym_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `term_id` int(10) unsigned NOT NULL,
  `name` text COLLATE utf8_general_ci NOT NULL,
  `type` enum('EXACT','BROAD','NARROW','RELATED') DEFAULT NULL,
  `dbxref` varchar(500) DEFAULT NULL,
  PRIMARY KEY (`synonym_id`),
  UNIQUE KEY `synonym_term_idx` (`term_id`,`synonym_id`),
  KEY `synonym_name_idx` (`name`(100))
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;

CREATE TABLE `alt_id` (
  `alt_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `term_id` int(10) unsigned NOT NULL,
  `accession` varchar(64) NOT NULL,
  PRIMARY KEY (`alt_id`),
  UNIQUE KEY `term_alt_idx` (`term_id`,`alt_id`),
  KEY `ix_alt_id_accession` (`accession`)
) ENGINE=MyISAM;

CREATE TABLE `relation_type` (
  `relation_type_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(64) NOT NULL,
  PRIMARY KEY (`relation_type_id`),
  UNIQUE KEY `name` (`name`)
) ENGINE=MyISAM;


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
) ENGINE=MyISAM;


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
) ENGINE=MyISAM;

-- There are additional tables in the released databases called
-- "aux_XX_YY_map".  These are created by the "add_subset_maps.pl"
-- scripts.  Please see the README document for further information.

-- $Id$
