-- MySQL dump 10.13  Distrib 5.1.61, for redhat-linux-gnu (x86_64)
--
-- Host: mysql-eg-devel-1.ebi.ac.uk    Database: homo_sapiens_core_test_db
-- ------------------------------------------------------
-- Server version	5.1.49-log

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `object_xref`
--

DROP TABLE IF EXISTS `object_xref`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `object_xref` (
  `object_xref_id` int(11) NOT NULL AUTO_INCREMENT,
  `ensembl_id` int(10) unsigned NOT NULL DEFAULT '0',
  `ensembl_object_type` enum('RawContig','Transcript','Gene','Translation','regulatory_factor','regulatory_feature') COLLATE latin1_bin NOT NULL DEFAULT 'RawContig',
  `xref_id` int(10) unsigned NOT NULL,
  `linkage_annotation` varchar(255) COLLATE latin1_bin DEFAULT NULL,
  `analysis_id` smallint(5) unsigned NOT NULL,
  UNIQUE KEY `ensembl_object_type` (`ensembl_object_type`,`ensembl_id`,`xref_id`),
  KEY `oxref_idx` (`object_xref_id`,`xref_id`,`ensembl_object_type`,`ensembl_id`),
  KEY `xref_idx` (`xref_id`,`ensembl_object_type`)
) ENGINE=MyISAM AUTO_INCREMENT=253685 DEFAULT CHARSET=latin1 COLLATE=latin1_bin;
/*!40101 SET character_set_client = @saved_cs_client */;

/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2012-11-05 10:52:11
