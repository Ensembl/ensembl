-- MySQL dump 10.13  Distrib 5.5.49, for debian-linux-gnu (x86_64)
--
-- Host: localhost    Database: ens_test
-- ------------------------------------------------------
-- Server version	5.5.49-0ubuntu0.14.04.1

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `allele`
--

DROP TABLE IF EXISTS `allele`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `allele` (
  `allele_id` int(11) NOT NULL AUTO_INCREMENT,
  `variation_id` int(11) unsigned NOT NULL,
  `subsnp_id` int(11) unsigned DEFAULT NULL,
  `allele_code_id` int(11) unsigned NOT NULL,
  `population_id` int(11) unsigned DEFAULT NULL,
  `frequency` float unsigned DEFAULT NULL,
  `count` int(11) unsigned DEFAULT NULL,
  `frequency_submitter_handle` int(10) DEFAULT NULL,
  PRIMARY KEY (`allele_id`),
  KEY `variation_idx` (`variation_id`),
  KEY `subsnp_idx` (`subsnp_id`),
  KEY `population_idx` (`population_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `allele`
--

LOCK TABLES `allele` WRITE;
/*!40000 ALTER TABLE `allele` DISABLE KEYS */;
/*!40000 ALTER TABLE `allele` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `allele_code`
--

DROP TABLE IF EXISTS `allele_code`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `allele_code` (
  `allele_code_id` int(11) NOT NULL AUTO_INCREMENT,
  `allele` varchar(60000) DEFAULT NULL,
  PRIMARY KEY (`allele_code_id`),
  UNIQUE KEY `allele_idx` (`allele`(1000))
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `allele_code`
--

LOCK TABLES `allele_code` WRITE;
/*!40000 ALTER TABLE `allele_code` DISABLE KEYS */;
/*!40000 ALTER TABLE `allele_code` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `associate_study`
--

DROP TABLE IF EXISTS `associate_study`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `associate_study` (
  `study1_id` int(10) unsigned NOT NULL,
  `study2_id` int(10) unsigned NOT NULL,
  PRIMARY KEY (`study1_id`,`study2_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `associate_study`
--

LOCK TABLES `associate_study` WRITE;
/*!40000 ALTER TABLE `associate_study` DISABLE KEYS */;
/*!40000 ALTER TABLE `associate_study` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `attrib`
--

DROP TABLE IF EXISTS `attrib`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `attrib` (
  `attrib_id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `attrib_type_id` smallint(5) unsigned NOT NULL DEFAULT '0',
  `value` text NOT NULL,
  PRIMARY KEY (`attrib_id`),
  UNIQUE KEY `type_val_idx` (`attrib_type_id`,`value`(80))
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `attrib`
--

LOCK TABLES `attrib` WRITE;
/*!40000 ALTER TABLE `attrib` DISABLE KEYS */;
/*!40000 ALTER TABLE `attrib` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `attrib_set`
--

DROP TABLE IF EXISTS `attrib_set`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `attrib_set` (
  `attrib_set_id` int(11) unsigned NOT NULL DEFAULT '0',
  `attrib_id` int(11) unsigned NOT NULL DEFAULT '0',
  UNIQUE KEY `set_idx` (`attrib_set_id`,`attrib_id`),
  KEY `attrib_idx` (`attrib_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `attrib_set`
--

LOCK TABLES `attrib_set` WRITE;
/*!40000 ALTER TABLE `attrib_set` DISABLE KEYS */;
/*!40000 ALTER TABLE `attrib_set` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `attrib_type`
--

DROP TABLE IF EXISTS `attrib_type`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `attrib_type` (
  `attrib_type_id` smallint(5) unsigned NOT NULL DEFAULT '0',
  `code` varchar(20) NOT NULL DEFAULT '',
  `name` varchar(255) NOT NULL DEFAULT '',
  `description` text,
  PRIMARY KEY (`attrib_type_id`),
  UNIQUE KEY `code_idx` (`code`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `attrib_type`
--

LOCK TABLES `attrib_type` WRITE;
/*!40000 ALTER TABLE `attrib_type` DISABLE KEYS */;
/*!40000 ALTER TABLE `attrib_type` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `compressed_genotype_region`
--

DROP TABLE IF EXISTS `compressed_genotype_region`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `compressed_genotype_region` (
  `sample_id` int(10) unsigned NOT NULL,
  `seq_region_id` int(10) unsigned NOT NULL,
  `seq_region_start` int(11) NOT NULL,
  `seq_region_end` int(11) NOT NULL,
  `seq_region_strand` tinyint(4) NOT NULL,
  `genotypes` blob,
  KEY `pos_idx` (`seq_region_id`,`seq_region_start`),
  KEY `sample_idx` (`sample_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `compressed_genotype_region`
--

LOCK TABLES `compressed_genotype_region` WRITE;
/*!40000 ALTER TABLE `compressed_genotype_region` DISABLE KEYS */;
/*!40000 ALTER TABLE `compressed_genotype_region` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `compressed_genotype_var`
--

DROP TABLE IF EXISTS `compressed_genotype_var`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `compressed_genotype_var` (
  `variation_id` int(11) unsigned NOT NULL,
  `subsnp_id` int(11) unsigned DEFAULT NULL,
  `genotypes` blob,
  KEY `variation_idx` (`variation_id`),
  KEY `subsnp_idx` (`subsnp_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `compressed_genotype_var`
--

LOCK TABLES `compressed_genotype_var` WRITE;
/*!40000 ALTER TABLE `compressed_genotype_var` DISABLE KEYS */;
/*!40000 ALTER TABLE `compressed_genotype_var` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `coord_system`
--

DROP TABLE IF EXISTS `coord_system`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `coord_system` (
  `coord_system_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `species_id` int(10) unsigned NOT NULL DEFAULT '1',
  `name` varchar(40) NOT NULL,
  `version` varchar(255) DEFAULT NULL,
  `rank` int(11) NOT NULL,
  `attrib` set('default_version','sequence_level') DEFAULT NULL,
  PRIMARY KEY (`coord_system_id`),
  UNIQUE KEY `rank_idx` (`rank`,`species_id`),
  UNIQUE KEY `name_idx` (`name`,`version`,`species_id`),
  KEY `species_idx` (`species_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `coord_system`
--

LOCK TABLES `coord_system` WRITE;
/*!40000 ALTER TABLE `coord_system` DISABLE KEYS */;
/*!40000 ALTER TABLE `coord_system` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `display_group`
--

DROP TABLE IF EXISTS `display_group`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `display_group` (
  `display_group_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `display_priority` int(10) unsigned NOT NULL,
  `display_name` varchar(255) NOT NULL,
  PRIMARY KEY (`display_group_id`),
  UNIQUE KEY `display_name` (`display_name`),
  UNIQUE KEY `display_priority` (`display_priority`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `display_group`
--

LOCK TABLES `display_group` WRITE;
/*!40000 ALTER TABLE `display_group` DISABLE KEYS */;
/*!40000 ALTER TABLE `display_group` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `failed_allele`
--

DROP TABLE IF EXISTS `failed_allele`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `failed_allele` (
  `failed_allele_id` int(11) NOT NULL AUTO_INCREMENT,
  `allele_id` int(10) unsigned NOT NULL,
  `failed_description_id` int(10) unsigned NOT NULL,
  PRIMARY KEY (`failed_allele_id`),
  UNIQUE KEY `allele_idx` (`allele_id`,`failed_description_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `failed_allele`
--

LOCK TABLES `failed_allele` WRITE;
/*!40000 ALTER TABLE `failed_allele` DISABLE KEYS */;
/*!40000 ALTER TABLE `failed_allele` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `failed_description`
--

DROP TABLE IF EXISTS `failed_description`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `failed_description` (
  `failed_description_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `description` text NOT NULL,
  PRIMARY KEY (`failed_description_id`)
) ENGINE=MyISAM AUTO_INCREMENT=21 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `failed_description`
--

LOCK TABLES `failed_description` WRITE;
/*!40000 ALTER TABLE `failed_description` DISABLE KEYS */;
INSERT INTO `failed_description` VALUES (1,'Variant maps to more than 3 different locations'),(2,'None of the variant alleles match the reference allele'),(3,'Variant has more than 3 different alleles'),(4,'Loci with no observed variant alleles in dbSNP'),(5,'Variant does not map to the genome'),(6,'Variant has no genotypes'),(7,'Genotype frequencies do not add up to 1'),(8,'Variant has no associated sequence'),(9,'Variant submission has been withdrawn by the 1000 genomes project due to high false positive rate'),(11,'Additional submitted allele data from dbSNP does not agree with the dbSNP refSNP alleles'),(12,'Variant has more than 3 different submitted alleles'),(13,'Alleles contain non-nucleotide characters'),(14,'Alleles contain ambiguity codes'),(15,'Mapped position is not compatible with reported alleles'),(16,'Flagged as suspect by dbSNP'),(17,'Variant can not be re-mapped to the current assembly'),(18,'Supporting evidence can not be re-mapped to the current assembly'),(19,'Variant maps to more than one genomic location'),(20,'Variant at first base in sequence');
/*!40000 ALTER TABLE `failed_description` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `failed_structural_variation`
--

DROP TABLE IF EXISTS `failed_structural_variation`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `failed_structural_variation` (
  `failed_structural_variation_id` int(11) NOT NULL AUTO_INCREMENT,
  `structural_variation_id` int(10) unsigned NOT NULL,
  `failed_description_id` int(10) unsigned NOT NULL,
  PRIMARY KEY (`failed_structural_variation_id`),
  UNIQUE KEY `structural_variation_idx` (`structural_variation_id`,`failed_description_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `failed_structural_variation`
--

LOCK TABLES `failed_structural_variation` WRITE;
/*!40000 ALTER TABLE `failed_structural_variation` DISABLE KEYS */;
/*!40000 ALTER TABLE `failed_structural_variation` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `failed_variation`
--

DROP TABLE IF EXISTS `failed_variation`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `failed_variation` (
  `failed_variation_id` int(11) NOT NULL AUTO_INCREMENT,
  `variation_id` int(10) unsigned NOT NULL,
  `failed_description_id` int(10) unsigned NOT NULL,
  PRIMARY KEY (`failed_variation_id`),
  UNIQUE KEY `variation_idx` (`variation_id`,`failed_description_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `failed_variation`
--

LOCK TABLES `failed_variation` WRITE;
/*!40000 ALTER TABLE `failed_variation` DISABLE KEYS */;
/*!40000 ALTER TABLE `failed_variation` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `genotype_code`
--

DROP TABLE IF EXISTS `genotype_code`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `genotype_code` (
  `genotype_code_id` int(11) unsigned NOT NULL,
  `allele_code_id` int(11) unsigned NOT NULL,
  `haplotype_id` tinyint(2) unsigned NOT NULL,
  `phased` tinyint(2) unsigned DEFAULT NULL,
  KEY `genotype_code_id` (`genotype_code_id`),
  KEY `allele_code_id` (`allele_code_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `genotype_code`
--

LOCK TABLES `genotype_code` WRITE;
/*!40000 ALTER TABLE `genotype_code` DISABLE KEYS */;
/*!40000 ALTER TABLE `genotype_code` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `individual`
--

DROP TABLE IF EXISTS `individual`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `individual` (
  `individual_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(255) DEFAULT NULL,
  `description` text,
  `gender` enum('Male','Female','Unknown') NOT NULL DEFAULT 'Unknown',
  `father_individual_id` int(10) unsigned DEFAULT NULL,
  `mother_individual_id` int(10) unsigned DEFAULT NULL,
  `individual_type_id` int(10) unsigned NOT NULL DEFAULT '0',
  PRIMARY KEY (`individual_id`),
  KEY `father_individual_idx` (`father_individual_id`),
  KEY `mother_individual_idx` (`mother_individual_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `individual`
--

LOCK TABLES `individual` WRITE;
/*!40000 ALTER TABLE `individual` DISABLE KEYS */;
/*!40000 ALTER TABLE `individual` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `individual_synonym`
--

DROP TABLE IF EXISTS `individual_synonym`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `individual_synonym` (
  `synonym_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `individual_id` int(10) unsigned NOT NULL,
  `source_id` int(10) unsigned NOT NULL,
  `name` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`synonym_id`),
  KEY `individual_idx` (`individual_id`),
  KEY `name` (`name`,`source_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `individual_synonym`
--

LOCK TABLES `individual_synonym` WRITE;
/*!40000 ALTER TABLE `individual_synonym` DISABLE KEYS */;
/*!40000 ALTER TABLE `individual_synonym` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `individual_type`
--

DROP TABLE IF EXISTS `individual_type`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `individual_type` (
  `individual_type_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(255) NOT NULL,
  `description` text,
  PRIMARY KEY (`individual_type_id`)
) ENGINE=MyISAM AUTO_INCREMENT=5 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `individual_type`
--

LOCK TABLES `individual_type` WRITE;
/*!40000 ALTER TABLE `individual_type` DISABLE KEYS */;
INSERT INTO `individual_type` VALUES (1,'fully_inbred','multiple organisms have the same genome sequence'),(2,'partly_inbred','single organisms have reduced genome variability due to human intervention'),(3,'outbred','a single organism which breeds freely'),(4,'mutant','a single or multiple organisms with the same genome sequence that have a natural or experimentally induced mutation');
/*!40000 ALTER TABLE `individual_type` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `meta`
--

DROP TABLE IF EXISTS `meta`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `meta` (
  `meta_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `species_id` int(10) unsigned DEFAULT '1',
  `meta_key` varchar(40) NOT NULL,
  `meta_value` varchar(255) NOT NULL,
  PRIMARY KEY (`meta_id`),
  UNIQUE KEY `species_key_value_idx` (`species_id`,`meta_key`,`meta_value`),
  KEY `species_value_idx` (`species_id`,`meta_value`)
) ENGINE=MyISAM AUTO_INCREMENT=10 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `meta`
--

LOCK TABLES `meta` WRITE;
/*!40000 ALTER TABLE `meta` DISABLE KEYS */;
INSERT INTO `meta` VALUES (1,NULL,'schema_type','variation'),(2,NULL,'schema_version','84'),(3,NULL,'patch','patch_84_85_a.sql|schema version'),(4,NULL,'patch','patch_84_85_b.sql|create sample_synonym'),(5,NULL,'patch','patch_84_85_c.sql|drop column moltype from variation_synonym'),(6,NULL,'patch','patch_84_85_d.sql|Making attrib_id auto_increment'),(7,NULL,'patch','patch_84_85_e.sql|drop the table tagged_variation_feature'),(8,NULL,'patch','patch_84_85_f.sql|add phenotype_ontology_accession'),(9,NULL,'patch','patch_84_85_g.sql|allow the column description to store more text in the source table');
/*!40000 ALTER TABLE `meta` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `meta_coord`
--

DROP TABLE IF EXISTS `meta_coord`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `meta_coord` (
  `table_name` varchar(40) NOT NULL,
  `coord_system_id` int(10) unsigned NOT NULL,
  `max_length` int(11) DEFAULT NULL,
  UNIQUE KEY `table_name` (`table_name`,`coord_system_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `meta_coord`
--

LOCK TABLES `meta_coord` WRITE;
/*!40000 ALTER TABLE `meta_coord` DISABLE KEYS */;
/*!40000 ALTER TABLE `meta_coord` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `motif_feature_variation`
--

DROP TABLE IF EXISTS `motif_feature_variation`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `motif_feature_variation` (
  `motif_feature_variation_id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `variation_feature_id` int(11) unsigned NOT NULL,
  `feature_stable_id` varchar(128) DEFAULT NULL,
  `motif_feature_id` int(11) unsigned NOT NULL,
  `allele_string` text,
  `somatic` tinyint(1) NOT NULL DEFAULT '0',
  `consequence_types` set('TF_binding_site_variant','TFBS_ablation','TFBS_fusion','TFBS_amplification','TFBS_translocation') DEFAULT NULL,
  `motif_name` varchar(60) DEFAULT NULL,
  `motif_start` int(11) unsigned DEFAULT NULL,
  `motif_end` int(11) unsigned DEFAULT NULL,
  `motif_score_delta` float DEFAULT NULL,
  `in_informative_position` tinyint(1) NOT NULL DEFAULT '0',
  PRIMARY KEY (`motif_feature_variation_id`),
  KEY `variation_feature_idx` (`variation_feature_id`),
  KEY `consequence_type_idx` (`consequence_types`),
  KEY `somatic_feature_idx` (`feature_stable_id`,`somatic`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `motif_feature_variation`
--

LOCK TABLES `motif_feature_variation` WRITE;
/*!40000 ALTER TABLE `motif_feature_variation` DISABLE KEYS */;
/*!40000 ALTER TABLE `motif_feature_variation` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `phenotype`
--

DROP TABLE IF EXISTS `phenotype`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `phenotype` (
  `phenotype_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `stable_id` varchar(255) DEFAULT NULL,
  `name` varchar(50) DEFAULT NULL,
  `description` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`phenotype_id`),
  UNIQUE KEY `desc_idx` (`description`),
  KEY `name_idx` (`name`),
  KEY `stable_idx` (`stable_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `phenotype`
--

LOCK TABLES `phenotype` WRITE;
/*!40000 ALTER TABLE `phenotype` DISABLE KEYS */;
/*!40000 ALTER TABLE `phenotype` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `phenotype_feature`
--

DROP TABLE IF EXISTS `phenotype_feature`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `phenotype_feature` (
  `phenotype_feature_id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `phenotype_id` int(11) unsigned DEFAULT NULL,
  `source_id` int(11) unsigned DEFAULT NULL,
  `study_id` int(11) unsigned DEFAULT NULL,
  `type` enum('Gene','Variation','StructuralVariation','SupportingStructuralVariation','QTL','RegulatoryFeature') DEFAULT NULL,
  `object_id` varchar(255) DEFAULT NULL,
  `is_significant` tinyint(1) unsigned DEFAULT '1',
  `seq_region_id` int(11) unsigned DEFAULT NULL,
  `seq_region_start` int(11) unsigned DEFAULT NULL,
  `seq_region_end` int(11) unsigned DEFAULT NULL,
  `seq_region_strand` tinyint(4) DEFAULT NULL,
  PRIMARY KEY (`phenotype_feature_id`),
  KEY `phenotype_idx` (`phenotype_id`),
  KEY `object_idx` (`object_id`,`type`),
  KEY `type_idx` (`type`),
  KEY `pos_idx` (`seq_region_id`,`seq_region_start`,`seq_region_end`),
  KEY `source_idx` (`source_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `phenotype_feature`
--

LOCK TABLES `phenotype_feature` WRITE;
/*!40000 ALTER TABLE `phenotype_feature` DISABLE KEYS */;
/*!40000 ALTER TABLE `phenotype_feature` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `phenotype_feature_attrib`
--

DROP TABLE IF EXISTS `phenotype_feature_attrib`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `phenotype_feature_attrib` (
  `phenotype_feature_id` int(11) unsigned NOT NULL,
  `attrib_type_id` int(11) DEFAULT NULL,
  `value` varchar(255) DEFAULT NULL,
  KEY `phenotype_feature_idx` (`phenotype_feature_id`),
  KEY `type_value_idx` (`attrib_type_id`,`value`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `phenotype_feature_attrib`
--

LOCK TABLES `phenotype_feature_attrib` WRITE;
/*!40000 ALTER TABLE `phenotype_feature_attrib` DISABLE KEYS */;
/*!40000 ALTER TABLE `phenotype_feature_attrib` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `phenotype_ontology_accession`
--

DROP TABLE IF EXISTS `phenotype_ontology_accession`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `phenotype_ontology_accession` (
  `phenotype_id` int(11) unsigned NOT NULL,
  `accession` varchar(255) NOT NULL,
  `linked_by_attrib` set('437','438','439','440','441','442','443','444') DEFAULT NULL,
  PRIMARY KEY (`phenotype_id`,`accession`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `phenotype_ontology_accession`
--

LOCK TABLES `phenotype_ontology_accession` WRITE;
/*!40000 ALTER TABLE `phenotype_ontology_accession` DISABLE KEYS */;
/*!40000 ALTER TABLE `phenotype_ontology_accession` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `population`
--

DROP TABLE IF EXISTS `population`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `population` (
  `population_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(255) DEFAULT NULL,
  `size` int(10) DEFAULT NULL,
  `description` text,
  `collection` tinyint(1) DEFAULT '0',
  `freqs_from_gts` tinyint(1) DEFAULT NULL,
  `display` enum('LD','MARTDISPLAYABLE','UNDISPLAYABLE') DEFAULT 'UNDISPLAYABLE',
  `display_group_id` tinyint(1) DEFAULT NULL,
  PRIMARY KEY (`population_id`),
  KEY `name_idx` (`name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `population`
--

LOCK TABLES `population` WRITE;
/*!40000 ALTER TABLE `population` DISABLE KEYS */;
/*!40000 ALTER TABLE `population` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `population_genotype`
--

DROP TABLE IF EXISTS `population_genotype`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `population_genotype` (
  `population_genotype_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `variation_id` int(11) unsigned NOT NULL,
  `subsnp_id` int(11) unsigned DEFAULT NULL,
  `genotype_code_id` int(11) DEFAULT NULL,
  `frequency` float DEFAULT NULL,
  `population_id` int(10) unsigned DEFAULT NULL,
  `count` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`population_genotype_id`),
  KEY `population_idx` (`population_id`),
  KEY `variation_idx` (`variation_id`),
  KEY `subsnp_idx` (`subsnp_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `population_genotype`
--

LOCK TABLES `population_genotype` WRITE;
/*!40000 ALTER TABLE `population_genotype` DISABLE KEYS */;
/*!40000 ALTER TABLE `population_genotype` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `population_structure`
--

DROP TABLE IF EXISTS `population_structure`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `population_structure` (
  `super_population_id` int(10) unsigned NOT NULL,
  `sub_population_id` int(10) unsigned NOT NULL,
  UNIQUE KEY `super_population_idx` (`super_population_id`,`sub_population_id`),
  KEY `sub_population_idx` (`sub_population_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `population_structure`
--

LOCK TABLES `population_structure` WRITE;
/*!40000 ALTER TABLE `population_structure` DISABLE KEYS */;
/*!40000 ALTER TABLE `population_structure` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `population_synonym`
--

DROP TABLE IF EXISTS `population_synonym`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `population_synonym` (
  `synonym_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `population_id` int(10) unsigned NOT NULL,
  `source_id` int(10) unsigned NOT NULL,
  `name` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`synonym_id`),
  KEY `population_idx` (`population_id`),
  KEY `name` (`name`,`source_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `population_synonym`
--

LOCK TABLES `population_synonym` WRITE;
/*!40000 ALTER TABLE `population_synonym` DISABLE KEYS */;
/*!40000 ALTER TABLE `population_synonym` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `protein_function_predictions`
--

DROP TABLE IF EXISTS `protein_function_predictions`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `protein_function_predictions` (
  `translation_md5_id` int(11) unsigned NOT NULL,
  `analysis_attrib_id` int(11) unsigned NOT NULL,
  `prediction_matrix` mediumblob,
  PRIMARY KEY (`translation_md5_id`,`analysis_attrib_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `protein_function_predictions`
--

LOCK TABLES `protein_function_predictions` WRITE;
/*!40000 ALTER TABLE `protein_function_predictions` DISABLE KEYS */;
/*!40000 ALTER TABLE `protein_function_predictions` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `protein_function_predictions_attrib`
--

DROP TABLE IF EXISTS `protein_function_predictions_attrib`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `protein_function_predictions_attrib` (
  `translation_md5_id` int(11) unsigned NOT NULL,
  `analysis_attrib_id` int(11) unsigned NOT NULL,
  `attrib_type_id` int(11) unsigned NOT NULL,
  `position_values` blob,
  PRIMARY KEY (`translation_md5_id`,`analysis_attrib_id`,`attrib_type_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `protein_function_predictions_attrib`
--

LOCK TABLES `protein_function_predictions_attrib` WRITE;
/*!40000 ALTER TABLE `protein_function_predictions_attrib` DISABLE KEYS */;
/*!40000 ALTER TABLE `protein_function_predictions_attrib` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `publication`
--

DROP TABLE IF EXISTS `publication`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `publication` (
  `publication_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `title` varchar(255) DEFAULT NULL,
  `authors` varchar(255) DEFAULT NULL,
  `pmid` int(10) DEFAULT NULL,
  `pmcid` varchar(255) DEFAULT NULL,
  `year` int(10) unsigned DEFAULT NULL,
  `doi` varchar(50) DEFAULT NULL,
  `ucsc_id` varchar(50) DEFAULT NULL,
  PRIMARY KEY (`publication_id`),
  KEY `pmid_idx` (`pmid`),
  KEY `doi_idx` (`doi`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `publication`
--

LOCK TABLES `publication` WRITE;
/*!40000 ALTER TABLE `publication` DISABLE KEYS */;
/*!40000 ALTER TABLE `publication` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `read_coverage`
--

DROP TABLE IF EXISTS `read_coverage`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `read_coverage` (
  `seq_region_id` int(10) unsigned NOT NULL,
  `seq_region_start` int(11) NOT NULL,
  `seq_region_end` int(11) NOT NULL,
  `level` tinyint(4) NOT NULL,
  `sample_id` int(10) unsigned NOT NULL,
  KEY `seq_region_idx` (`seq_region_id`,`seq_region_start`),
  KEY `sample_idx` (`sample_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `read_coverage`
--

LOCK TABLES `read_coverage` WRITE;
/*!40000 ALTER TABLE `read_coverage` DISABLE KEYS */;
/*!40000 ALTER TABLE `read_coverage` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `regulatory_feature_variation`
--

DROP TABLE IF EXISTS `regulatory_feature_variation`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `regulatory_feature_variation` (
  `regulatory_feature_variation_id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `variation_feature_id` int(11) unsigned NOT NULL,
  `feature_stable_id` varchar(128) DEFAULT NULL,
  `feature_type` text,
  `allele_string` text,
  `somatic` tinyint(1) NOT NULL DEFAULT '0',
  `consequence_types` set('regulatory_region_variant','regulatory_region_ablation','regulatory_region_fusion','regulatory_region_amplification','regulatory_region_translocation') DEFAULT NULL,
  PRIMARY KEY (`regulatory_feature_variation_id`),
  KEY `variation_feature_idx` (`variation_feature_id`),
  KEY `consequence_type_idx` (`consequence_types`),
  KEY `somatic_feature_idx` (`feature_stable_id`,`somatic`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `regulatory_feature_variation`
--

LOCK TABLES `regulatory_feature_variation` WRITE;
/*!40000 ALTER TABLE `regulatory_feature_variation` DISABLE KEYS */;
/*!40000 ALTER TABLE `regulatory_feature_variation` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `sample`
--

DROP TABLE IF EXISTS `sample`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `sample` (
  `sample_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `individual_id` int(10) unsigned NOT NULL,
  `name` varchar(255) DEFAULT NULL,
  `description` text,
  `study_id` int(10) unsigned DEFAULT NULL,
  `display` enum('REFERENCE','DEFAULT','DISPLAYABLE','UNDISPLAYABLE','LD','MARTDISPLAYABLE') DEFAULT 'UNDISPLAYABLE',
  `has_coverage` tinyint(1) unsigned NOT NULL DEFAULT '0',
  `variation_set_id` set('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31','32','33','34','35','36','37','38','39','40','41','42','43','44','45','46','47','48','49','50','51','52','53','54','55','56','57','58','59','60','61','62','63','64') DEFAULT NULL,
  PRIMARY KEY (`sample_id`),
  KEY `individual_idx` (`individual_id`),
  KEY `study_idx` (`study_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `sample`
--

LOCK TABLES `sample` WRITE;
/*!40000 ALTER TABLE `sample` DISABLE KEYS */;
/*!40000 ALTER TABLE `sample` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `sample_genotype_multiple_bp`
--

DROP TABLE IF EXISTS `sample_genotype_multiple_bp`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `sample_genotype_multiple_bp` (
  `variation_id` int(10) unsigned NOT NULL,
  `subsnp_id` int(15) unsigned DEFAULT NULL,
  `allele_1` varchar(25000) DEFAULT NULL,
  `allele_2` varchar(25000) DEFAULT NULL,
  `sample_id` int(10) unsigned DEFAULT NULL,
  KEY `variation_idx` (`variation_id`),
  KEY `subsnp_idx` (`subsnp_id`),
  KEY `sample_idx` (`sample_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `sample_genotype_multiple_bp`
--

LOCK TABLES `sample_genotype_multiple_bp` WRITE;
/*!40000 ALTER TABLE `sample_genotype_multiple_bp` DISABLE KEYS */;
/*!40000 ALTER TABLE `sample_genotype_multiple_bp` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `sample_population`
--

DROP TABLE IF EXISTS `sample_population`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `sample_population` (
  `sample_id` int(10) unsigned NOT NULL,
  `population_id` int(10) unsigned NOT NULL,
  KEY `sample_idx` (`sample_id`),
  KEY `population_idx` (`population_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `sample_population`
--

LOCK TABLES `sample_population` WRITE;
/*!40000 ALTER TABLE `sample_population` DISABLE KEYS */;
/*!40000 ALTER TABLE `sample_population` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `sample_synonym`
--

DROP TABLE IF EXISTS `sample_synonym`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `sample_synonym` (
  `synonym_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `sample_id` int(10) unsigned NOT NULL,
  `source_id` int(10) unsigned NOT NULL,
  `name` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`synonym_id`),
  KEY `sample_idx` (`sample_id`),
  KEY `name` (`name`,`source_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `sample_synonym`
--

LOCK TABLES `sample_synonym` WRITE;
/*!40000 ALTER TABLE `sample_synonym` DISABLE KEYS */;
/*!40000 ALTER TABLE `sample_synonym` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `seq_region`
--

DROP TABLE IF EXISTS `seq_region`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `seq_region` (
  `seq_region_id` int(10) unsigned NOT NULL,
  `name` varchar(40) NOT NULL,
  `coord_system_id` int(10) unsigned NOT NULL,
  PRIMARY KEY (`seq_region_id`),
  UNIQUE KEY `name_cs_idx` (`name`,`coord_system_id`),
  KEY `cs_idx` (`coord_system_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `seq_region`
--

LOCK TABLES `seq_region` WRITE;
/*!40000 ALTER TABLE `seq_region` DISABLE KEYS */;
/*!40000 ALTER TABLE `seq_region` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `source`
--

DROP TABLE IF EXISTS `source`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `source` (
  `source_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(24) NOT NULL,
  `version` int(11) DEFAULT NULL,
  `description` varchar(400) DEFAULT NULL,
  `url` varchar(255) DEFAULT NULL,
  `type` enum('chip','lsdb') DEFAULT NULL,
  `somatic_status` enum('germline','somatic','mixed') DEFAULT 'germline',
  `data_types` set('variation','variation_synonym','structural_variation','phenotype_feature','study') DEFAULT NULL,
  PRIMARY KEY (`source_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `source`
--

LOCK TABLES `source` WRITE;
/*!40000 ALTER TABLE `source` DISABLE KEYS */;
/*!40000 ALTER TABLE `source` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `strain_gtype_poly`
--

DROP TABLE IF EXISTS `strain_gtype_poly`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `strain_gtype_poly` (
  `variation_id` int(10) unsigned NOT NULL,
  `sample_name` varchar(100) DEFAULT NULL,
  PRIMARY KEY (`variation_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `strain_gtype_poly`
--

LOCK TABLES `strain_gtype_poly` WRITE;
/*!40000 ALTER TABLE `strain_gtype_poly` DISABLE KEYS */;
/*!40000 ALTER TABLE `strain_gtype_poly` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `structural_variation`
--

DROP TABLE IF EXISTS `structural_variation`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `structural_variation` (
  `structural_variation_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `variation_name` varchar(255) DEFAULT NULL,
  `alias` varchar(255) DEFAULT NULL,
  `source_id` int(10) unsigned NOT NULL,
  `study_id` int(10) unsigned DEFAULT NULL,
  `class_attrib_id` int(10) unsigned NOT NULL DEFAULT '0',
  `clinical_significance` set('uncertain significance','not provided','benign','likely benign','likely pathogenic','pathogenic','drug response','histocompatibility','other','confers sensitivity','risk factor','association','protective') DEFAULT NULL,
  `validation_status` enum('validated','not validated','high quality') DEFAULT NULL,
  `is_evidence` tinyint(4) DEFAULT '0',
  `somatic` tinyint(1) NOT NULL DEFAULT '0',
  `copy_number` tinyint(2) DEFAULT NULL,
  PRIMARY KEY (`structural_variation_id`),
  UNIQUE KEY `variation_name` (`variation_name`),
  KEY `source_idx` (`source_id`),
  KEY `study_idx` (`study_id`),
  KEY `attrib_idx` (`class_attrib_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `structural_variation`
--

LOCK TABLES `structural_variation` WRITE;
/*!40000 ALTER TABLE `structural_variation` DISABLE KEYS */;
/*!40000 ALTER TABLE `structural_variation` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `structural_variation_association`
--

DROP TABLE IF EXISTS `structural_variation_association`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `structural_variation_association` (
  `structural_variation_id` int(10) unsigned NOT NULL,
  `supporting_structural_variation_id` int(10) unsigned NOT NULL,
  PRIMARY KEY (`structural_variation_id`,`supporting_structural_variation_id`),
  KEY `structural_variation_idx` (`structural_variation_id`),
  KEY `supporting_structural_variation_idx` (`supporting_structural_variation_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `structural_variation_association`
--

LOCK TABLES `structural_variation_association` WRITE;
/*!40000 ALTER TABLE `structural_variation_association` DISABLE KEYS */;
/*!40000 ALTER TABLE `structural_variation_association` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `structural_variation_feature`
--

DROP TABLE IF EXISTS `structural_variation_feature`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `structural_variation_feature` (
  `structural_variation_feature_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `seq_region_id` int(10) unsigned NOT NULL,
  `outer_start` int(11) DEFAULT NULL,
  `seq_region_start` int(11) NOT NULL,
  `inner_start` int(11) DEFAULT NULL,
  `inner_end` int(11) DEFAULT NULL,
  `seq_region_end` int(11) NOT NULL,
  `outer_end` int(11) DEFAULT NULL,
  `seq_region_strand` tinyint(4) NOT NULL,
  `structural_variation_id` int(10) unsigned NOT NULL,
  `variation_name` varchar(255) DEFAULT NULL,
  `source_id` int(10) unsigned NOT NULL,
  `study_id` int(10) unsigned DEFAULT NULL,
  `class_attrib_id` int(10) unsigned NOT NULL DEFAULT '0',
  `allele_string` longtext,
  `is_evidence` tinyint(1) NOT NULL DEFAULT '0',
  `somatic` tinyint(1) NOT NULL DEFAULT '0',
  `breakpoint_order` tinyint(4) DEFAULT NULL,
  `length` int(10) DEFAULT NULL,
  `variation_set_id` set('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31','32','33','34','35','36','37','38','39','40','41','42','43','44','45','46','47','48','49','50','51','52','53','54','55','56','57','58','59','60','61','62','63','64') NOT NULL DEFAULT '',
  PRIMARY KEY (`structural_variation_feature_id`),
  KEY `pos_idx` (`seq_region_id`,`seq_region_start`,`seq_region_end`),
  KEY `structural_variation_idx` (`structural_variation_id`),
  KEY `source_idx` (`source_id`),
  KEY `study_idx` (`study_id`),
  KEY `attrib_idx` (`class_attrib_id`),
  KEY `variation_set_idx` (`variation_set_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `structural_variation_feature`
--

LOCK TABLES `structural_variation_feature` WRITE;
/*!40000 ALTER TABLE `structural_variation_feature` DISABLE KEYS */;
/*!40000 ALTER TABLE `structural_variation_feature` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `structural_variation_sample`
--

DROP TABLE IF EXISTS `structural_variation_sample`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `structural_variation_sample` (
  `structural_variation_sample_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `structural_variation_id` int(10) unsigned NOT NULL,
  `sample_id` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`structural_variation_sample_id`),
  KEY `structural_variation_idx` (`structural_variation_id`),
  KEY `sample_idx` (`sample_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `structural_variation_sample`
--

LOCK TABLES `structural_variation_sample` WRITE;
/*!40000 ALTER TABLE `structural_variation_sample` DISABLE KEYS */;
/*!40000 ALTER TABLE `structural_variation_sample` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `study`
--

DROP TABLE IF EXISTS `study`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `study` (
  `study_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `source_id` int(10) unsigned NOT NULL,
  `name` varchar(255) DEFAULT NULL,
  `description` text,
  `url` varchar(255) DEFAULT NULL,
  `external_reference` varchar(255) DEFAULT NULL,
  `study_type` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`study_id`),
  KEY `source_idx` (`source_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `study`
--

LOCK TABLES `study` WRITE;
/*!40000 ALTER TABLE `study` DISABLE KEYS */;
/*!40000 ALTER TABLE `study` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `submitter_handle`
--

DROP TABLE IF EXISTS `submitter_handle`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `submitter_handle` (
  `handle_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `handle` varchar(25) DEFAULT NULL,
  PRIMARY KEY (`handle_id`),
  UNIQUE KEY `handle` (`handle`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `submitter_handle`
--

LOCK TABLES `submitter_handle` WRITE;
/*!40000 ALTER TABLE `submitter_handle` DISABLE KEYS */;
/*!40000 ALTER TABLE `submitter_handle` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `subsnp_handle`
--

DROP TABLE IF EXISTS `subsnp_handle`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `subsnp_handle` (
  `subsnp_id` int(11) unsigned NOT NULL,
  `handle` varchar(20) DEFAULT NULL,
  PRIMARY KEY (`subsnp_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `subsnp_handle`
--

LOCK TABLES `subsnp_handle` WRITE;
/*!40000 ALTER TABLE `subsnp_handle` DISABLE KEYS */;
/*!40000 ALTER TABLE `subsnp_handle` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `tmp_sample_genotype_single_bp`
--

DROP TABLE IF EXISTS `tmp_sample_genotype_single_bp`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `tmp_sample_genotype_single_bp` (
  `variation_id` int(10) NOT NULL,
  `subsnp_id` int(15) unsigned DEFAULT NULL,
  `allele_1` char(1) DEFAULT NULL,
  `allele_2` char(1) DEFAULT NULL,
  `sample_id` int(10) unsigned NOT NULL,
  KEY `variation_idx` (`variation_id`),
  KEY `subsnp_idx` (`subsnp_id`),
  KEY `sample_idx` (`sample_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=100000000;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `tmp_sample_genotype_single_bp`
--

LOCK TABLES `tmp_sample_genotype_single_bp` WRITE;
/*!40000 ALTER TABLE `tmp_sample_genotype_single_bp` DISABLE KEYS */;
/*!40000 ALTER TABLE `tmp_sample_genotype_single_bp` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `transcript_variation`
--

DROP TABLE IF EXISTS `transcript_variation`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `transcript_variation` (
  `transcript_variation_id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `variation_feature_id` int(11) unsigned NOT NULL,
  `feature_stable_id` varchar(128) DEFAULT NULL,
  `allele_string` text,
  `somatic` tinyint(1) NOT NULL DEFAULT '0',
  `consequence_types` set('splice_acceptor_variant','splice_donor_variant','stop_lost','coding_sequence_variant','missense_variant','stop_gained','synonymous_variant','frameshift_variant','non_coding_transcript_variant','non_coding_transcript_exon_variant','mature_miRNA_variant','NMD_transcript_variant','5_prime_UTR_variant','3_prime_UTR_variant','incomplete_terminal_codon_variant','intron_variant','splice_region_variant','downstream_gene_variant','upstream_gene_variant','start_lost','stop_retained_variant','inframe_insertion','inframe_deletion','transcript_ablation','transcript_fusion','transcript_amplification','transcript_translocation','TFBS_ablation','TFBS_fusion','TFBS_amplification','TFBS_translocation','regulatory_region_ablation','regulatory_region_fusion','regulatory_region_amplification','regulatory_region_translocation','feature_elongation','feature_truncation','protein_altering_variant') DEFAULT NULL,
  `cds_start` int(11) unsigned DEFAULT NULL,
  `cds_end` int(11) unsigned DEFAULT NULL,
  `cdna_start` int(11) unsigned DEFAULT NULL,
  `cdna_end` int(11) unsigned DEFAULT NULL,
  `translation_start` int(11) unsigned DEFAULT NULL,
  `translation_end` int(11) unsigned DEFAULT NULL,
  `distance_to_transcript` int(11) unsigned DEFAULT NULL,
  `codon_allele_string` text,
  `pep_allele_string` text,
  `hgvs_genomic` text,
  `hgvs_transcript` text,
  `hgvs_protein` text,
  `polyphen_prediction` enum('unknown','benign','possibly damaging','probably damaging') DEFAULT NULL,
  `polyphen_score` float DEFAULT NULL,
  `sift_prediction` enum('tolerated','deleterious','tolerated - low confidence','deleterious - low confidence') DEFAULT NULL,
  `sift_score` float DEFAULT NULL,
  `display` int(1) DEFAULT '1',
  PRIMARY KEY (`transcript_variation_id`),
  KEY `variation_feature_idx` (`variation_feature_id`),
  KEY `consequence_type_idx` (`consequence_types`),
  KEY `somatic_feature_idx` (`feature_stable_id`,`somatic`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `transcript_variation`
--

LOCK TABLES `transcript_variation` WRITE;
/*!40000 ALTER TABLE `transcript_variation` DISABLE KEYS */;
/*!40000 ALTER TABLE `transcript_variation` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `translation_md5`
--

DROP TABLE IF EXISTS `translation_md5`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `translation_md5` (
  `translation_md5_id` int(11) NOT NULL AUTO_INCREMENT,
  `translation_md5` char(32) NOT NULL,
  PRIMARY KEY (`translation_md5_id`),
  UNIQUE KEY `md5_idx` (`translation_md5`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `translation_md5`
--

LOCK TABLES `translation_md5` WRITE;
/*!40000 ALTER TABLE `translation_md5` DISABLE KEYS */;
/*!40000 ALTER TABLE `translation_md5` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `variation`
--

DROP TABLE IF EXISTS `variation`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `variation` (
  `variation_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `source_id` int(10) unsigned NOT NULL,
  `name` varchar(255) DEFAULT NULL,
  `ancestral_allele` varchar(255) DEFAULT NULL,
  `flipped` tinyint(1) unsigned DEFAULT NULL,
  `class_attrib_id` int(10) unsigned DEFAULT '0',
  `somatic` tinyint(1) NOT NULL DEFAULT '0',
  `minor_allele` varchar(50) DEFAULT NULL,
  `minor_allele_freq` float DEFAULT NULL,
  `minor_allele_count` int(10) unsigned DEFAULT NULL,
  `clinical_significance` set('uncertain significance','not provided','benign','likely benign','likely pathogenic','pathogenic','drug response','histocompatibility','other','confers sensitivity','risk factor','association','protective') DEFAULT NULL,
  `evidence_attribs` set('367','368','369','370','371','372','418','421') DEFAULT NULL,
  `display` int(1) DEFAULT '1',
  PRIMARY KEY (`variation_id`),
  UNIQUE KEY `name` (`name`),
  KEY `source_idx` (`source_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `variation`
--

LOCK TABLES `variation` WRITE;
/*!40000 ALTER TABLE `variation` DISABLE KEYS */;
/*!40000 ALTER TABLE `variation` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `variation_attrib`
--

DROP TABLE IF EXISTS `variation_attrib`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `variation_attrib` (
  `variation_id` int(11) unsigned NOT NULL,
  `attrib_id` int(11) DEFAULT NULL,
  `value` varchar(255) DEFAULT NULL,
  KEY `variation_idx` (`variation_id`),
  KEY `attrib_value_idx` (`attrib_id`,`value`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `variation_attrib`
--

LOCK TABLES `variation_attrib` WRITE;
/*!40000 ALTER TABLE `variation_attrib` DISABLE KEYS */;
/*!40000 ALTER TABLE `variation_attrib` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `variation_citation`
--

DROP TABLE IF EXISTS `variation_citation`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `variation_citation` (
  `variation_id` int(10) unsigned NOT NULL,
  `publication_id` int(10) unsigned NOT NULL,
  PRIMARY KEY (`variation_id`,`publication_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `variation_citation`
--

LOCK TABLES `variation_citation` WRITE;
/*!40000 ALTER TABLE `variation_citation` DISABLE KEYS */;
/*!40000 ALTER TABLE `variation_citation` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `variation_feature`
--

DROP TABLE IF EXISTS `variation_feature`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `variation_feature` (
  `variation_feature_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `seq_region_id` int(10) unsigned NOT NULL,
  `seq_region_start` int(11) NOT NULL,
  `seq_region_end` int(11) NOT NULL,
  `seq_region_strand` tinyint(4) NOT NULL,
  `variation_id` int(10) unsigned NOT NULL,
  `allele_string` varchar(50000) DEFAULT NULL,
  `variation_name` varchar(255) DEFAULT NULL,
  `map_weight` int(11) NOT NULL,
  `flags` set('genotyped') DEFAULT NULL,
  `source_id` int(10) unsigned NOT NULL,
  `consequence_types` set('intergenic_variant','splice_acceptor_variant','splice_donor_variant','stop_lost','coding_sequence_variant','missense_variant','stop_gained','synonymous_variant','frameshift_variant','non_coding_transcript_variant','non_coding_transcript_exon_variant','mature_miRNA_variant','NMD_transcript_variant','5_prime_UTR_variant','3_prime_UTR_variant','incomplete_terminal_codon_variant','intron_variant','splice_region_variant','downstream_gene_variant','upstream_gene_variant','start_lost','stop_retained_variant','inframe_insertion','inframe_deletion','transcript_ablation','transcript_fusion','transcript_amplification','transcript_translocation','TFBS_ablation','TFBS_fusion','TFBS_amplification','TFBS_translocation','regulatory_region_ablation','regulatory_region_fusion','regulatory_region_amplification','regulatory_region_translocation','feature_elongation','feature_truncation','regulatory_region_variant','TF_binding_site_variant','protein_altering_variant') NOT NULL DEFAULT 'intergenic_variant',
  `variation_set_id` set('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31','32','33','34','35','36','37','38','39','40','41','42','43','44','45','46','47','48','49','50','51','52','53','54','55','56','57','58','59','60','61','62','63','64') NOT NULL DEFAULT '',
  `class_attrib_id` int(10) unsigned DEFAULT '0',
  `somatic` tinyint(1) NOT NULL DEFAULT '0',
  `minor_allele` varchar(50) DEFAULT NULL,
  `minor_allele_freq` float DEFAULT NULL,
  `minor_allele_count` int(10) unsigned DEFAULT NULL,
  `alignment_quality` double DEFAULT NULL,
  `evidence_attribs` set('367','368','369','370','371','372','418','421') DEFAULT NULL,
  `clinical_significance` set('uncertain significance','not provided','benign','likely benign','likely pathogenic','pathogenic','drug response','histocompatibility','other','confers sensitivity','risk factor','association','protective') DEFAULT NULL,
  `display` int(1) DEFAULT '1',
  PRIMARY KEY (`variation_feature_id`),
  KEY `pos_idx` (`seq_region_id`,`seq_region_start`,`seq_region_end`),
  KEY `variation_idx` (`variation_id`),
  KEY `variation_set_idx` (`variation_set_id`),
  KEY `consequence_type_idx` (`consequence_types`),
  KEY `source_idx` (`source_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `variation_feature`
--

LOCK TABLES `variation_feature` WRITE;
/*!40000 ALTER TABLE `variation_feature` DISABLE KEYS */;
/*!40000 ALTER TABLE `variation_feature` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `variation_genename`
--

DROP TABLE IF EXISTS `variation_genename`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `variation_genename` (
  `variation_id` int(10) unsigned NOT NULL,
  `gene_name` varchar(255) NOT NULL,
  PRIMARY KEY (`variation_id`,`gene_name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `variation_genename`
--

LOCK TABLES `variation_genename` WRITE;
/*!40000 ALTER TABLE `variation_genename` DISABLE KEYS */;
/*!40000 ALTER TABLE `variation_genename` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `variation_hgvs`
--

DROP TABLE IF EXISTS `variation_hgvs`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `variation_hgvs` (
  `variation_id` int(10) unsigned NOT NULL,
  `hgvs_name` varchar(255) NOT NULL,
  PRIMARY KEY (`variation_id`,`hgvs_name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `variation_hgvs`
--

LOCK TABLES `variation_hgvs` WRITE;
/*!40000 ALTER TABLE `variation_hgvs` DISABLE KEYS */;
/*!40000 ALTER TABLE `variation_hgvs` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `variation_set`
--

DROP TABLE IF EXISTS `variation_set`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `variation_set` (
  `variation_set_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(255) DEFAULT NULL,
  `description` text,
  `short_name_attrib_id` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`variation_set_id`),
  KEY `name_idx` (`name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `variation_set`
--

LOCK TABLES `variation_set` WRITE;
/*!40000 ALTER TABLE `variation_set` DISABLE KEYS */;
/*!40000 ALTER TABLE `variation_set` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `variation_set_structural_variation`
--

DROP TABLE IF EXISTS `variation_set_structural_variation`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `variation_set_structural_variation` (
  `structural_variation_id` int(10) unsigned NOT NULL,
  `variation_set_id` int(10) unsigned NOT NULL,
  PRIMARY KEY (`structural_variation_id`,`variation_set_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `variation_set_structural_variation`
--

LOCK TABLES `variation_set_structural_variation` WRITE;
/*!40000 ALTER TABLE `variation_set_structural_variation` DISABLE KEYS */;
/*!40000 ALTER TABLE `variation_set_structural_variation` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `variation_set_structure`
--

DROP TABLE IF EXISTS `variation_set_structure`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `variation_set_structure` (
  `variation_set_super` int(10) unsigned NOT NULL,
  `variation_set_sub` int(10) unsigned NOT NULL,
  PRIMARY KEY (`variation_set_super`,`variation_set_sub`),
  KEY `sub_idx` (`variation_set_sub`,`variation_set_super`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `variation_set_structure`
--

LOCK TABLES `variation_set_structure` WRITE;
/*!40000 ALTER TABLE `variation_set_structure` DISABLE KEYS */;
/*!40000 ALTER TABLE `variation_set_structure` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `variation_set_variation`
--

DROP TABLE IF EXISTS `variation_set_variation`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `variation_set_variation` (
  `variation_id` int(10) unsigned NOT NULL,
  `variation_set_id` int(10) unsigned NOT NULL,
  PRIMARY KEY (`variation_id`,`variation_set_id`),
  KEY `variation_set_idx` (`variation_set_id`,`variation_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `variation_set_variation`
--

LOCK TABLES `variation_set_variation` WRITE;
/*!40000 ALTER TABLE `variation_set_variation` DISABLE KEYS */;
/*!40000 ALTER TABLE `variation_set_variation` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `variation_synonym`
--

DROP TABLE IF EXISTS `variation_synonym`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `variation_synonym` (
  `variation_synonym_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `variation_id` int(10) unsigned NOT NULL,
  `subsnp_id` int(15) unsigned DEFAULT NULL,
  `source_id` int(10) unsigned NOT NULL,
  `name` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`variation_synonym_id`),
  UNIQUE KEY `name` (`name`,`source_id`),
  KEY `variation_idx` (`variation_id`),
  KEY `subsnp_idx` (`subsnp_id`),
  KEY `source_idx` (`source_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `variation_synonym`
--

LOCK TABLES `variation_synonym` WRITE;
/*!40000 ALTER TABLE `variation_synonym` DISABLE KEYS */;
/*!40000 ALTER TABLE `variation_synonym` ENABLE KEYS */;
UNLOCK TABLES;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2016-06-02 13:39:28
