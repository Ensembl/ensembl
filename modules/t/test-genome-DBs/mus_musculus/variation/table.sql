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
) ENGINE=MyISAM  ;

CREATE TABLE `allele_code` (
  `allele_code_id` int(11) NOT NULL AUTO_INCREMENT,
  `allele` varchar(100) DEFAULT NULL,
  PRIMARY KEY (`allele_code_id`),
  UNIQUE KEY `allele_idx` (`allele`)
) ENGINE=MyISAM  ;

CREATE TABLE `associate_study` (
  `study1_id` int(10) unsigned NOT NULL,
  `study2_id` int(10) unsigned NOT NULL,
  PRIMARY KEY (`study1_id`,`study2_id`)
) ENGINE=MyISAM ;

CREATE TABLE `attrib` (
  `attrib_id` int(11) unsigned NOT NULL DEFAULT '0',
  `attrib_type_id` smallint(5) unsigned NOT NULL DEFAULT '0',
  `value` text NOT NULL,
  PRIMARY KEY (`attrib_id`),
  UNIQUE KEY `type_val_idx` (`attrib_type_id`,`value`(80))
) ENGINE=MyISAM ;

CREATE TABLE `attrib_set` (
  `attrib_set_id` int(11) unsigned NOT NULL DEFAULT '0',
  `attrib_id` int(11) unsigned NOT NULL DEFAULT '0',
  UNIQUE KEY `set_idx` (`attrib_set_id`,`attrib_id`),
  KEY `attrib_idx` (`attrib_id`)
) ENGINE=MyISAM ;

CREATE TABLE `attrib_type` (
  `attrib_type_id` smallint(5) unsigned NOT NULL DEFAULT '0',
  `code` varchar(20) NOT NULL DEFAULT '',
  `name` varchar(255) NOT NULL DEFAULT '',
  `description` text,
  PRIMARY KEY (`attrib_type_id`),
  UNIQUE KEY `code_idx` (`code`)
) ENGINE=MyISAM ;

CREATE TABLE `compressed_genotype_region` (
  `sample_id` int(10) unsigned NOT NULL,
  `seq_region_id` int(10) unsigned NOT NULL,
  `seq_region_start` int(11) NOT NULL,
  `seq_region_end` int(11) NOT NULL,
  `seq_region_strand` tinyint(4) NOT NULL,
  `genotypes` blob,
  KEY `pos_idx` (`seq_region_id`,`seq_region_start`),
  KEY `sample_idx` (`sample_id`)
) ENGINE=MyISAM ;

CREATE TABLE `compressed_genotype_var` (
  `variation_id` int(11) unsigned NOT NULL,
  `subsnp_id` int(11) unsigned DEFAULT NULL,
  `genotypes` blob,
  KEY `variation_idx` (`variation_id`),
  KEY `subsnp_idx` (`subsnp_id`)
) ENGINE=MyISAM ;

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
) ENGINE=MyISAM ;

CREATE TABLE `display_group` (
  `display_group_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `display_priority` int(10) unsigned NOT NULL,
  `display_name` varchar(255) NOT NULL,
  PRIMARY KEY (`display_group_id`),
  UNIQUE KEY `display_name` (`display_name`),
  UNIQUE KEY `display_priority` (`display_priority`)
) ENGINE=MyISAM  ;

CREATE TABLE `failed_allele` (
  `failed_allele_id` int(11) NOT NULL AUTO_INCREMENT,
  `allele_id` int(10) unsigned NOT NULL,
  `failed_description_id` int(10) unsigned NOT NULL,
  PRIMARY KEY (`failed_allele_id`),
  UNIQUE KEY `allele_idx` (`allele_id`,`failed_description_id`)
) ENGINE=MyISAM  ;

CREATE TABLE `failed_description` (
  `failed_description_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `description` text NOT NULL,
  PRIMARY KEY (`failed_description_id`)
) ENGINE=MyISAM  ;

CREATE TABLE `failed_structural_variation` (
  `failed_structural_variation_id` int(11) NOT NULL AUTO_INCREMENT,
  `structural_variation_id` int(10) unsigned NOT NULL,
  `failed_description_id` int(10) unsigned NOT NULL,
  PRIMARY KEY (`failed_structural_variation_id`),
  UNIQUE KEY `structural_variation_idx` (`structural_variation_id`,`failed_description_id`)
) ENGINE=MyISAM  ;

CREATE TABLE `failed_variation` (
  `failed_variation_id` int(11) NOT NULL AUTO_INCREMENT,
  `variation_id` int(10) unsigned NOT NULL,
  `failed_description_id` int(10) unsigned NOT NULL,
  PRIMARY KEY (`failed_variation_id`),
  UNIQUE KEY `variation_idx` (`variation_id`,`failed_description_id`)
) ENGINE=MyISAM  ;

CREATE TABLE `genotype_code` (
  `genotype_code_id` int(11) unsigned NOT NULL,
  `allele_code_id` int(11) unsigned NOT NULL,
  `haplotype_id` tinyint(2) unsigned NOT NULL,
  `phased` tinyint(2) unsigned DEFAULT NULL,
  KEY `genotype_code_id` (`genotype_code_id`),
  KEY `allele_code_id` (`allele_code_id`)
) ENGINE=MyISAM ;

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
) ENGINE=MyISAM  ;

CREATE TABLE `individual_synonym` (
  `synonym_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `individual_id` int(10) unsigned NOT NULL,
  `source_id` int(10) unsigned NOT NULL,
  `name` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`synonym_id`),
  KEY `individual_idx` (`individual_id`),
  KEY `name` (`name`,`source_id`)
) ENGINE=MyISAM  ;

CREATE TABLE `individual_type` (
  `individual_type_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(255) NOT NULL,
  `description` text,
  PRIMARY KEY (`individual_type_id`)
) ENGINE=MyISAM  ;

CREATE TABLE `meta` (
  `meta_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `species_id` int(10) unsigned DEFAULT '1',
  `meta_key` varchar(40) NOT NULL,
  `meta_value` varchar(255) NOT NULL,
  PRIMARY KEY (`meta_id`),
  UNIQUE KEY `species_key_value_idx` (`species_id`,`meta_key`,`meta_value`),
  KEY `species_value_idx` (`species_id`,`meta_value`)
) ENGINE=MyISAM  ;

CREATE TABLE `meta_coord` (
  `table_name` varchar(40) NOT NULL,
  `coord_system_id` int(10) unsigned NOT NULL,
  `max_length` int(11) DEFAULT NULL,
  UNIQUE KEY `table_name` (`table_name`,`coord_system_id`)
) ENGINE=MyISAM ;

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
) ENGINE=MyISAM  ;

CREATE TABLE `phenotype` (
  `phenotype_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `stable_id` varchar(255) DEFAULT NULL,
  `name` varchar(50) DEFAULT NULL,
  `description` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`phenotype_id`),
  UNIQUE KEY `desc_idx` (`description`),
  KEY `name_idx` (`name`),
  KEY `stable_idx` (`stable_id`)
) ENGINE=MyISAM  ;

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
) ENGINE=MyISAM  ;

CREATE TABLE `phenotype_feature_attrib` (
  `phenotype_feature_id` int(11) unsigned NOT NULL,
  `attrib_type_id` int(11) DEFAULT NULL,
  `value` varchar(255) DEFAULT NULL,
  KEY `phenotype_feature_idx` (`phenotype_feature_id`),
  KEY `type_value_idx` (`attrib_type_id`,`value`)
) ENGINE=MyISAM ;

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
) ENGINE=MyISAM  ;

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
) ENGINE=MyISAM  ;

CREATE TABLE `population_structure` (
  `super_population_id` int(10) unsigned NOT NULL,
  `sub_population_id` int(10) unsigned NOT NULL,
  UNIQUE KEY `super_population_idx` (`super_population_id`,`sub_population_id`),
  KEY `sub_population_idx` (`sub_population_id`)
) ENGINE=MyISAM ;

CREATE TABLE `population_synonym` (
  `synonym_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `population_id` int(10) unsigned NOT NULL,
  `source_id` int(10) unsigned NOT NULL,
  `name` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`synonym_id`),
  KEY `population_idx` (`population_id`),
  KEY `name` (`name`,`source_id`)
) ENGINE=MyISAM  ;

CREATE TABLE `protein_function_predictions` (
  `translation_md5_id` int(11) unsigned NOT NULL,
  `analysis_attrib_id` int(11) unsigned NOT NULL,
  `prediction_matrix` mediumblob,
  PRIMARY KEY (`translation_md5_id`,`analysis_attrib_id`)
) ENGINE=MyISAM ;

CREATE TABLE `protein_function_predictions_attrib` (
  `translation_md5_id` int(11) unsigned NOT NULL,
  `analysis_attrib_id` int(11) unsigned NOT NULL,
  `attrib_type_id` int(11) unsigned NOT NULL,
  `position_values` blob,
  PRIMARY KEY (`translation_md5_id`,`analysis_attrib_id`,`attrib_type_id`)
) ENGINE=MyISAM ;

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
) ENGINE=MyISAM  ;

CREATE TABLE `read_coverage` (
  `seq_region_id` int(10) unsigned NOT NULL,
  `seq_region_start` int(11) NOT NULL,
  `seq_region_end` int(11) NOT NULL,
  `level` tinyint(4) NOT NULL,
  `sample_id` int(10) unsigned NOT NULL,
  KEY `seq_region_idx` (`seq_region_id`,`seq_region_start`),
  KEY `sample_idx` (`sample_id`)
) ENGINE=MyISAM ;

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
) ENGINE=MyISAM  ;

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
) ENGINE=MyISAM  ;

CREATE TABLE `sample_genotype_multiple_bp` (
  `variation_id` int(10) unsigned NOT NULL,
  `subsnp_id` int(15) unsigned DEFAULT NULL,
  `allele_1` varchar(25000) DEFAULT NULL,
  `allele_2` varchar(25000) DEFAULT NULL,
  `sample_id` int(10) unsigned DEFAULT NULL,
  KEY `variation_idx` (`variation_id`),
  KEY `subsnp_idx` (`subsnp_id`),
  KEY `sample_idx` (`sample_id`)
) ENGINE=MyISAM ;

CREATE TABLE `sample_population` (
  `sample_id` int(10) unsigned NOT NULL,
  `population_id` int(10) unsigned NOT NULL,
  KEY `population_idx` (`population_id`),
  KEY `sample_idx` (`sample_id`)
) ENGINE=MyISAM ;

CREATE TABLE `seq_region` (
  `seq_region_id` int(10) unsigned NOT NULL,
  `name` varchar(40) NOT NULL,
  `coord_system_id` int(10) unsigned NOT NULL,
  PRIMARY KEY (`seq_region_id`),
  UNIQUE KEY `name_cs_idx` (`name`,`coord_system_id`),
  KEY `cs_idx` (`coord_system_id`)
) ENGINE=MyISAM ;

CREATE TABLE `source` (
  `source_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(24) NOT NULL,
  `version` int(11) DEFAULT NULL,
  `description` varchar(255) DEFAULT NULL,
  `url` varchar(255) DEFAULT NULL,
  `type` enum('chip','lsdb') DEFAULT NULL,
  `somatic_status` enum('germline','somatic','mixed') DEFAULT 'germline',
  `data_types` set('variation','variation_synonym','structural_variation','phenotype_feature','study') DEFAULT NULL,
  PRIMARY KEY (`source_id`)
) ENGINE=MyISAM  ;

CREATE TABLE `strain_gtype_poly` (
  `variation_id` int(10) unsigned NOT NULL,
  `CAST_EI` varchar(100) DEFAULT NULL,
  `NOD_LTJ` varchar(100) DEFAULT NULL,
  `C57BL_10J` varchar(100) DEFAULT NULL,
  `A_J` varchar(100) DEFAULT NULL,
  `129_NBT` varchar(100) DEFAULT NULL,
  `129X1_SVJ` varchar(100) DEFAULT NULL,
  `MOLF_EIJ` varchar(100) DEFAULT NULL,
  `PL_J` varchar(100) DEFAULT NULL,
  `129S5` varchar(100) DEFAULT NULL,
  `129S4_SVJAE` varchar(100) DEFAULT NULL,
  `FVB` varchar(100) DEFAULT NULL,
  `strain_4` varchar(100) DEFAULT NULL,
  `LCAD_B6NTAC_DEFICIENT` varchar(100) DEFAULT NULL,
  `CAST_EIJ` varchar(100) DEFAULT NULL,
  `LG_J` varchar(100) DEFAULT NULL,
  `129P3_J` varchar(100) DEFAULT NULL,
  `CASA_RKJ` varchar(100) DEFAULT NULL,
  `SPRET_EIJ` varchar(100) DEFAULT NULL,
  `SEG_PAS` varchar(100) DEFAULT NULL,
  `MSM_MS` varchar(100) DEFAULT NULL,
  `WSB_EIJ` varchar(100) DEFAULT NULL,
  `PERA_EIJ` varchar(100) DEFAULT NULL,
  `BTBR_T__TF_J` varchar(100) DEFAULT NULL,
  `NON_LTJ` varchar(100) DEFAULT NULL,
  `CZECHII_EI` varchar(100) DEFAULT NULL,
  `AJ` varchar(100) DEFAULT NULL,
  `I_LNJ` varchar(100) DEFAULT NULL,
  `DBA_2J` varchar(100) DEFAULT NULL,
  `LPJ` varchar(100) DEFAULT NULL,
  `PERA_EI` varchar(100) DEFAULT NULL,
  `PWK_PHJ` varchar(100) DEFAULT NULL,
  `NZW_LACJ` varchar(100) DEFAULT NULL,
  `MRL_MPJ` varchar(100) DEFAULT NULL,
  `BLG2_MS` varchar(100) DEFAULT NULL,
  `DBA_1J` varchar(100) DEFAULT NULL,
  `CHD_MS` varchar(100) DEFAULT NULL,
  `BALB_CBY` varchar(100) DEFAULT NULL,
  `SM_J` varchar(100) DEFAULT NULL,
  `NZB_BLNJ` varchar(100) DEFAULT NULL,
  `A_HE` varchar(100) DEFAULT NULL,
  `MOLF` varchar(100) DEFAULT NULL,
  `SEA_GNJ` varchar(100) DEFAULT NULL,
  `PERA` varchar(100) DEFAULT NULL,
  `strain_3` varchar(100) DEFAULT NULL,
  `DBA_2` varchar(100) DEFAULT NULL,
  `SPRETUS_EIJ` varchar(100) DEFAULT NULL,
  `CZECHII_EIJ` varchar(100) DEFAULT NULL,
  `C3H_HE` varchar(100) DEFAULT NULL,
  `RF_J` varchar(100) DEFAULT NULL,
  `PWD_PH` varchar(100) DEFAULT NULL,
  `SJL_J` varchar(100) DEFAULT NULL,
  `NZB_BINJ` varchar(100) DEFAULT NULL,
  `CSS_NBT` varchar(100) DEFAULT NULL,
  `129P2` varchar(100) DEFAULT NULL,
  `SPRET_EI` varchar(100) DEFAULT NULL,
  `C58_J` varchar(100) DEFAULT NULL,
  `ST_BJ` varchar(100) DEFAULT NULL,
  `BALBC` varchar(100) DEFAULT NULL,
  `PWD_PHJ` varchar(100) DEFAULT NULL,
  `SAMR1` varchar(100) DEFAULT NULL,
  `SAMP8` varchar(100) DEFAULT NULL,
  `NZB_BLN` varchar(100) DEFAULT NULL,
  `SRET_EI` varchar(100) DEFAULT NULL,
  `A` varchar(100) DEFAULT NULL,
  `HMI_MS` varchar(100) DEFAULT NULL,
  `B10_D2` varchar(100) DEFAULT NULL,
  `BALB_CUCSD` varchar(100) DEFAULT NULL,
  `NZW_LAC` varchar(100) DEFAULT NULL,
  `BALB_C` varchar(100) DEFAULT NULL,
  `MA_MYJ` varchar(100) DEFAULT NULL,
  `JF1_MS` varchar(100) DEFAULT NULL,
  `RIIIS_J` varchar(100) DEFAULT NULL,
  `SJL_NBOM` varchar(100) DEFAULT NULL,
  `MSM` varchar(100) DEFAULT NULL,
  `C57L_J` varchar(100) DEFAULT NULL,
  `129X1SVJ` varchar(100) DEFAULT NULL,
  `SMJ` varchar(100) DEFAULT NULL,
  `C57BR_CDJ` varchar(100) DEFAULT NULL,
  `BALB_CA` varchar(100) DEFAULT NULL,
  `NOD_MRKTAC` varchar(100) DEFAULT NULL,
  `NOD__SHILTJ` varchar(100) DEFAULT NULL,
  `BUB_BNJ` varchar(100) DEFAULT NULL,
  `NOD_DIL` varchar(100) DEFAULT NULL,
  `A_HEJ` varchar(100) DEFAULT NULL,
  `DDK_PAS` varchar(100) DEFAULT NULL,
  `C3H_HEJ` varchar(100) DEFAULT NULL,
  `KK_HLJ` varchar(100) DEFAULT NULL,
  `NOD_SHILTJ` varchar(100) DEFAULT NULL,
  `129_SV` varchar(100) DEFAULT NULL,
  `CE_J` varchar(100) DEFAULT NULL,
  `MRL_MP` varchar(100) DEFAULT NULL,
  `BFM_2MS` varchar(100) DEFAULT NULL,
  `AKR` varchar(100) DEFAULT NULL,
  `129S6SVEV` varchar(100) DEFAULT NULL,
  `strain_5` varchar(100) DEFAULT NULL,
  `BUB` varchar(100) DEFAULT NULL,
  `B10_D2_H2_OSNJ` varchar(100) DEFAULT NULL,
  `MAI_PAS` varchar(100) DEFAULT NULL,
  `NJL_MS` varchar(100) DEFAULT NULL,
  `ZALENDE_EIJ` varchar(100) DEFAULT NULL,
  `C57BLKS_J` varchar(100) DEFAULT NULL,
  `SPRETUS` varchar(100) DEFAULT NULL,
  `SAMP1` varchar(100) DEFAULT NULL,
  `SWN_MS` varchar(100) DEFAULT NULL,
  `129P3J` varchar(100) DEFAULT NULL,
  `129X1_SV` varchar(100) DEFAULT NULL,
  `MUS_MUSCULUS_C57BL_6J` varchar(100) DEFAULT NULL,
  `NZO_HILTJ` varchar(100) DEFAULT NULL,
  `MUS_SPRETUS_SPRET_EI` varchar(100) DEFAULT NULL,
  `TSOD` varchar(100) DEFAULT NULL,
  `strain_1` varchar(100) DEFAULT NULL,
  `PGN2_MS` varchar(100) DEFAULT NULL,
  `KJR_MS` varchar(100) DEFAULT NULL,
  `AVZ_MS` varchar(100) DEFAULT NULL,
  `B10_D2_HC0_H2D_H2_T18C_OSNJ` varchar(100) DEFAULT NULL,
  `LGJ` varchar(100) DEFAULT NULL,
  `ERO_129P4_J` varchar(100) DEFAULT NULL,
  `ERO_SPRET_Ei` varchar(100) DEFAULT NULL,
  `ERO_C57BL_6J` varchar(100) DEFAULT NULL,
  `PI_DPPM_INT_C3H_HeJ` varchar(100) DEFAULT NULL,
  `PI_DPPM_INT_SWR_J` varchar(100) DEFAULT NULL,
  `JAX_129S1_SVIMJ` varchar(100) DEFAULT NULL,
  `PI_DPPM_INT_BALB_cJ` varchar(100) DEFAULT NULL,
  `ERO_BALB_cByJ` varchar(100) DEFAULT NULL,
  `ERO_AKR` varchar(100) DEFAULT NULL,
  `SC_MOUSE_GENOMES_CBA_J` varchar(100) DEFAULT NULL,
  `SC_MOUSE_GENOMES_AKR_J` varchar(100) DEFAULT NULL,
  `ERO_AKR_J` varchar(100) DEFAULT NULL,
  `SC_MOUSE_GENOMES_LP_J` varchar(100) DEFAULT NULL,
  `ERO_CBA_J` varchar(100) DEFAULT NULL,
  `ERO_C3H_HeJ` varchar(100) DEFAULT NULL,
  `ERO_129S6_SvEvTac` varchar(100) DEFAULT NULL,
  `ERO_A_J` varchar(100) DEFAULT NULL,
  `FVB_Ntac` varchar(100) DEFAULT NULL,
  `C57BL_6ByJ` varchar(100) DEFAULT NULL,
  `wild_17` varchar(100) DEFAULT NULL,
  `C57BL_6lcr` varchar(100) DEFAULT NULL,
  `PL_J_3` varchar(100) DEFAULT NULL,
  `C57BL_6HanZtm` varchar(100) DEFAULT NULL,
  `wild_7` varchar(100) DEFAULT NULL,
  `C57BL_6OlaHsdF` varchar(100) DEFAULT NULL,
  `C57BL_6OlaHsDc` varchar(100) DEFAULT NULL,
  `C57BL_6OlaHsdG` varchar(100) DEFAULT NULL,
  `wild_19` varchar(100) DEFAULT NULL,
  `C57BL_6Crl` varchar(100) DEFAULT NULL,
  `C57BL_6NSim` varchar(100) DEFAULT NULL,
  `C57BL_6N_e84_` varchar(100) DEFAULT NULL,
  `C57BL_6NTacAce` varchar(100) DEFAULT NULL,
  `BALB_cByJ_2` varchar(100) DEFAULT NULL,
  `BALB_cByJ_4` varchar(100) DEFAULT NULL,
  `wild_13` varchar(100) DEFAULT NULL,
  `DBA_2J_3` varchar(100) DEFAULT NULL,
  `C57BL_6BomTac` varchar(100) DEFAULT NULL,
  `PL_J_1` varchar(100) DEFAULT NULL,
  `SJL_J_1` varchar(100) DEFAULT NULL,
  `C3H_HeJ_1` varchar(100) DEFAULT NULL,
  `DBA_2J_1` varchar(100) DEFAULT NULL,
  `C3H_HeJ_2` varchar(100) DEFAULT NULL,
  `AKR_J_1` varchar(100) DEFAULT NULL,
  `C57BL_6JEiJ` varchar(100) DEFAULT NULL,
  `DBA_2J_2` varchar(100) DEFAULT NULL,
  `LP_J_1` varchar(100) DEFAULT NULL,
  `wild_11` varchar(100) DEFAULT NULL,
  `A_J_1` varchar(100) DEFAULT NULL,
  `AKR_J_3` varchar(100) DEFAULT NULL,
  `C3HeB_FeJ_3` varchar(100) DEFAULT NULL,
  `C3HeB_FeJ_2` varchar(100) DEFAULT NULL,
  `129S2_3` varchar(100) DEFAULT NULL,
  `PL_J_2` varchar(100) DEFAULT NULL,
  `wild_8` varchar(100) DEFAULT NULL,
  `A_J_3` varchar(100) DEFAULT NULL,
  `wild_15` varchar(100) DEFAULT NULL,
  `wild_9` varchar(100) DEFAULT NULL,
  `LP_J_2` varchar(100) DEFAULT NULL,
  `wild_12` varchar(100) DEFAULT NULL,
  `C3H_HeJ_3` varchar(100) DEFAULT NULL,
  `129S2_2` varchar(100) DEFAULT NULL,
  `wild_1` varchar(100) DEFAULT NULL,
  `wild_21` varchar(100) DEFAULT NULL,
  `LP_J_3` varchar(100) DEFAULT NULL,
  `AKR_J_2` varchar(100) DEFAULT NULL,
  `BALB_cByJ_1` varchar(100) DEFAULT NULL,
  `wild_10` varchar(100) DEFAULT NULL,
  `wild_5` varchar(100) DEFAULT NULL,
  `wild_6` varchar(100) DEFAULT NULL,
  `wild_14` varchar(100) DEFAULT NULL,
  `SJL_J_3` varchar(100) DEFAULT NULL,
  `wild_16` varchar(100) DEFAULT NULL,
  `129S2_4` varchar(100) DEFAULT NULL,
  `wild_18` varchar(100) DEFAULT NULL,
  `wild_2` varchar(100) DEFAULT NULL,
  `C57BL_6J_3` varchar(100) DEFAULT NULL,
  `SJL_J_2` varchar(100) DEFAULT NULL,
  `A_J_2` varchar(100) DEFAULT NULL,
  `129S2_1` varchar(100) DEFAULT NULL,
  `BALB_cByJ_3` varchar(100) DEFAULT NULL,
  `wild_4` varchar(100) DEFAULT NULL,
  `C3HeB_FeJ_1` varchar(100) DEFAULT NULL,
  `Cast_Ei_3` varchar(100) DEFAULT NULL,
  `Cast_Ei_2` varchar(100) DEFAULT NULL,
  `wild_3` varchar(100) DEFAULT NULL,
  `C57BL_6J_1` varchar(100) DEFAULT NULL,
  `C57BL_6J_2` varchar(100) DEFAULT NULL,
  `Cast_Ei_1` varchar(100) DEFAULT NULL,
  `C57BL_6J_5` varchar(100) DEFAULT NULL,
  `Spret_Ei_1` varchar(100) DEFAULT NULL,
  `Spret_Ei_2` varchar(100) DEFAULT NULL,
  `Spret_Ei_3` varchar(100) DEFAULT NULL,
  `wild_20` varchar(100) DEFAULT NULL,
  `C57BL_6J_4` varchar(100) DEFAULT NULL,
  `DBA2J` varchar(100) DEFAULT NULL,
  `CBAJ` varchar(100) DEFAULT NULL,
  `BALBcJ` varchar(100) DEFAULT NULL,
  `NZOHlLtJ` varchar(100) DEFAULT NULL,
  `C57BL6NJ` varchar(100) DEFAULT NULL,
  `129P2OlaHsd` varchar(100) DEFAULT NULL,
  `SPRETEiJ` varchar(100) DEFAULT NULL,
  `AKRJ` varchar(100) DEFAULT NULL,
  `PWKPhJ` varchar(100) DEFAULT NULL,
  `WSBEiJ` varchar(100) DEFAULT NULL,
  `129S5SvEvBrd` varchar(100) DEFAULT NULL,
  `NODShiLtJ` varchar(100) DEFAULT NULL,
  `129S1SvImJ` varchar(100) DEFAULT NULL,
  `CASTEiJ` varchar(100) DEFAULT NULL,
  `C3HHeJ` varchar(100) DEFAULT NULL,
  `A18_P11` varchar(100) DEFAULT NULL,
  `B86_P6` varchar(100) DEFAULT NULL,
  `A390_P7` varchar(100) DEFAULT NULL,
  `B91_P9` varchar(100) DEFAULT NULL,
  `B96_P8` varchar(100) DEFAULT NULL,
  `A292_P8` varchar(100) DEFAULT NULL,
  `B76_P13` varchar(100) DEFAULT NULL,
  `A202_P6` varchar(100) DEFAULT NULL,
  `231_Trig` varchar(100) DEFAULT NULL,
  `A382_P9` varchar(100) DEFAULT NULL,
  `A394_P5` varchar(100) DEFAULT NULL,
  `Mouse_C57BL_6` varchar(100) DEFAULT NULL,
  `Mouse_C57BL_6_x_SJL` varchar(100) DEFAULT NULL,
  `129S2` varchar(100) DEFAULT NULL,
  `C57BL_6N_Fam70b_tm1a_EUCOMM_Hmgu__H` varchar(100) DEFAULT NULL,
  `C57BL_6Brd_C57BL_6JIco_Tyr_c_Brd_` varchar(100) DEFAULT NULL,
  `B6NTac_B6N_Hint1_tm1a_EUCOMM_Wtsi__H` varchar(100) DEFAULT NULL,
  `B6NTac_B6N_Arhgef4_tm1a_KOMP_Wtsi__Wtsi` varchar(100) DEFAULT NULL,
  `C57BL_6N_Mmp12_tm1a_EUCOMM_Hmgu__H` varchar(100) DEFAULT NULL,
  `C57BL_6NTac_Ccdc101_tm1a_EUCOMM_Hmgu__Cnrm` varchar(100) DEFAULT NULL,
  `129S5_129P2` varchar(100) DEFAULT NULL,
  `C57BL_6Tac` varchar(100) DEFAULT NULL,
  `129S8_SvEv_Gpi1_c__NimrH` varchar(100) DEFAULT NULL,
  `Balb_c_129S2` varchar(100) DEFAULT NULL,
  `reference_C57BL_6J` varchar(100) DEFAULT NULL,
  PRIMARY KEY (`variation_id`)
) ENGINE=MyISAM ;

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
  KEY `study_idx` (`study_id`),
  KEY `attrib_idx` (`class_attrib_id`),
  KEY `source_idx` (`source_id`)
) ENGINE=MyISAM  ;

CREATE TABLE `structural_variation_association` (
  `structural_variation_id` int(10) unsigned NOT NULL,
  `supporting_structural_variation_id` int(10) unsigned NOT NULL,
  PRIMARY KEY (`structural_variation_id`,`supporting_structural_variation_id`),
  KEY `structural_variation_idx` (`structural_variation_id`),
  KEY `supporting_structural_variation_idx` (`supporting_structural_variation_id`)
) ENGINE=MyISAM ;

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
  `variation_set_id` set('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31','32','33','34','35','36','37','38','39','40','41','42','43','44','45','46','47','48','49','50','51','52','53','54','55','56','57','58','59','60','61','62','63','64') NOT NULL DEFAULT '',
  `somatic` tinyint(1) NOT NULL DEFAULT '0',
  `breakpoint_order` tinyint(4) DEFAULT NULL,
  `length` int(10) DEFAULT NULL,
  PRIMARY KEY (`structural_variation_feature_id`),
  KEY `pos_idx` (`seq_region_id`,`seq_region_start`,`seq_region_end`),
  KEY `structural_variation_idx` (`structural_variation_id`),
  KEY `attrib_idx` (`class_attrib_id`),
  KEY `source_idx` (`source_id`),
  KEY `variation_set_idx` (`variation_set_id`),
  KEY `study_idx` (`study_id`)
) ENGINE=MyISAM  ;

CREATE TABLE `structural_variation_sample` (
  `structural_variation_sample_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `structural_variation_id` int(10) unsigned NOT NULL,
  `sample_id` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`structural_variation_sample_id`),
  KEY `structural_variation_idx` (`structural_variation_id`),
  KEY `sample_idx` (`sample_id`)
) ENGINE=MyISAM  ;

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
) ENGINE=MyISAM  ;

CREATE TABLE `submitter_handle` (
  `handle_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `handle` varchar(25) DEFAULT NULL,
  PRIMARY KEY (`handle_id`),
  UNIQUE KEY `handle` (`handle`)
) ENGINE=MyISAM  ;

CREATE TABLE `subsnp_handle` (
  `subsnp_id` int(11) unsigned NOT NULL,
  `handle` varchar(20) DEFAULT NULL,
  PRIMARY KEY (`subsnp_id`)
) ENGINE=MyISAM ;

CREATE TABLE `subsnp_map` (
  `variation_id` int(11) unsigned NOT NULL,
  `subsnp_id` int(11) unsigned DEFAULT NULL,
  KEY `variation_idx` (`variation_id`)
) ENGINE=MyISAM ;

CREATE TABLE `tagged_variation_feature` (
  `variation_feature_id` int(10) unsigned NOT NULL,
  `tagged_variation_feature_id` int(10) unsigned DEFAULT NULL,
  `population_id` int(10) unsigned NOT NULL,
  KEY `tag_idx` (`variation_feature_id`),
  KEY `tagged_idx` (`tagged_variation_feature_id`),
  KEY `population_idx` (`population_id`)
) ENGINE=MyISAM ;

CREATE TABLE `tmp_sample_genotype_single_bp` (
  `variation_id` int(10) NOT NULL,
  `subsnp_id` int(15) unsigned DEFAULT NULL,
  `allele_1` char(1) DEFAULT NULL,
  `allele_2` char(1) DEFAULT NULL,
  `sample_id` int(10) unsigned NOT NULL,
  KEY `variation_idx` (`variation_id`),
  KEY `subsnp_idx` (`subsnp_id`),
  KEY `sample_idx` (`sample_id`)
) ENGINE=MyISAM  ;

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
  `sift_prediction` enum('tolerated','deleterious','tolerated - low confidence','deleterious - low confidence') DEFAULT NULL,
  `polyphen_score` float DEFAULT NULL,
  `sift_score` float DEFAULT NULL,
  `display` int(1) DEFAULT '1',
  PRIMARY KEY (`transcript_variation_id`),
  KEY `variation_feature_idx` (`variation_feature_id`),
  KEY `consequence_type_idx` (`consequence_types`),
  KEY `somatic_feature_idx` (`feature_stable_id`,`somatic`)
) ENGINE=MyISAM  ;

CREATE TABLE `translation_md5` (
  `translation_md5_id` int(11) NOT NULL AUTO_INCREMENT,
  `translation_md5` char(32) NOT NULL,
  PRIMARY KEY (`translation_md5_id`),
  UNIQUE KEY `md5_idx` (`translation_md5`)
) ENGINE=MyISAM  ;

CREATE TABLE `variation` (
  `variation_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `source_id` int(10) unsigned NOT NULL,
  `name` varchar(255) DEFAULT NULL,
  `validation_status` set('cluster','freq','submitter','doublehit','hapmap','1000Genome','failed','precious') DEFAULT NULL,
  `ancestral_allele` varchar(255) DEFAULT NULL,
  `flipped` tinyint(1) unsigned DEFAULT NULL,
  `class_attrib_id` int(10) unsigned DEFAULT '0',
  `somatic` tinyint(1) NOT NULL DEFAULT '0',
  `minor_allele` varchar(50) DEFAULT NULL,
  `minor_allele_freq` float DEFAULT NULL,
  `minor_allele_count` int(10) unsigned DEFAULT NULL,
  `clinical_significance` set('uncertain significance','not provided','benign','likely benign','likely pathogenic','pathogenic','drug response','histocompatibility','other','confers sensitivity','risk factor','association','protective') DEFAULT NULL,
  `evidence_attribs` set('367','368','369','370','371','372','418') DEFAULT NULL,
  `display` int(1) DEFAULT '1',
  PRIMARY KEY (`variation_id`),
  UNIQUE KEY `name` (`name`),
  KEY `source_idx` (`source_id`)
) ENGINE=MyISAM  ;

CREATE TABLE `variation_attrib` (
  `variation_id` int(11) unsigned NOT NULL,
  `attrib_id` int(11) DEFAULT NULL,
  `value` varchar(255) DEFAULT NULL,
  KEY `variation_idx` (`variation_id`),
  KEY `attrib_value_idx` (`attrib_id`,`value`)
) ENGINE=MyISAM ;

CREATE TABLE `variation_citation` (
  `variation_id` int(10) unsigned NOT NULL,
  `publication_id` int(10) unsigned NOT NULL,
  PRIMARY KEY (`variation_id`,`publication_id`)
) ENGINE=MyISAM ;

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
  `validation_status` set('cluster','freq','submitter','doublehit','hapmap','1000Genome','precious') DEFAULT NULL,
  `consequence_types` set('intergenic_variant','splice_acceptor_variant','splice_donor_variant','stop_lost','coding_sequence_variant','missense_variant','stop_gained','synonymous_variant','frameshift_variant','non_coding_transcript_variant','non_coding_transcript_exon_variant','mature_miRNA_variant','NMD_transcript_variant','5_prime_UTR_variant','3_prime_UTR_variant','incomplete_terminal_codon_variant','intron_variant','splice_region_variant','downstream_gene_variant','upstream_gene_variant','start_lost','stop_retained_variant','inframe_insertion','inframe_deletion','transcript_ablation','transcript_fusion','transcript_amplification','transcript_translocation','TFBS_ablation','TFBS_fusion','TFBS_amplification','TFBS_translocation','regulatory_region_ablation','regulatory_region_fusion','regulatory_region_amplification','regulatory_region_translocation','feature_elongation','feature_truncation','regulatory_region_variant','TF_binding_site_variant','protein_altering_variant') NOT NULL DEFAULT 'intergenic_variant',
  `variation_set_id` set('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31','32','33','34','35','36','37','38','39','40','41','42','43','44','45','46','47','48','49','50','51','52','53','54','55','56','57','58','59','60','61','62','63','64') NOT NULL DEFAULT '',
  `class_attrib_id` int(10) unsigned DEFAULT '0',
  `somatic` tinyint(1) NOT NULL DEFAULT '0',
  `minor_allele` varchar(50) DEFAULT NULL,
  `minor_allele_freq` float DEFAULT NULL,
  `minor_allele_count` int(10) unsigned DEFAULT NULL,
  `alignment_quality` double DEFAULT NULL,
  `evidence_attribs` set('367','368','369','370','371','372','418') DEFAULT NULL,
  `clinical_significance` set('uncertain significance','not provided','benign','likely benign','likely pathogenic','pathogenic','drug response','histocompatibility','other','confers sensitivity','risk factor','association','protective') DEFAULT NULL,
  `display` int(1) DEFAULT '1',
  PRIMARY KEY (`variation_feature_id`),
  KEY `pos_idx` (`seq_region_id`,`seq_region_start`,`seq_region_end`),
  KEY `variation_idx` (`variation_id`),
  KEY `variation_set_idx` (`variation_set_id`),
  KEY `consequence_type_idx` (`consequence_types`),
  KEY `source_idx` (`source_id`)
) ENGINE=MyISAM  ;

CREATE TABLE `variation_genename` (
  `variation_id` int(10) unsigned NOT NULL,
  `gene_name` varchar(255) NOT NULL,
  PRIMARY KEY (`variation_id`,`gene_name`)
) ENGINE=MyISAM ;

CREATE TABLE `variation_hgvs` (
  `variation_id` int(10) unsigned NOT NULL,
  `hgvs_name` varchar(255) NOT NULL,
  PRIMARY KEY (`variation_id`,`hgvs_name`)
) ENGINE=MyISAM ;

CREATE TABLE `variation_set` (
  `variation_set_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(255) DEFAULT NULL,
  `description` text,
  `short_name_attrib_id` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`variation_set_id`),
  KEY `name_idx` (`name`)
) ENGINE=MyISAM  ;

CREATE TABLE `variation_set_structural_variation` (
  `structural_variation_id` int(10) unsigned NOT NULL,
  `variation_set_id` int(10) unsigned NOT NULL,
  PRIMARY KEY (`structural_variation_id`,`variation_set_id`)
) ENGINE=MyISAM ;

CREATE TABLE `variation_set_structure` (
  `variation_set_super` int(10) unsigned NOT NULL,
  `variation_set_sub` int(10) unsigned NOT NULL,
  PRIMARY KEY (`variation_set_super`,`variation_set_sub`),
  KEY `sub_idx` (`variation_set_sub`,`variation_set_super`)
) ENGINE=MyISAM ;

CREATE TABLE `variation_set_variation` (
  `variation_id` int(10) unsigned NOT NULL,
  `variation_set_id` int(10) unsigned NOT NULL,
  PRIMARY KEY (`variation_id`,`variation_set_id`),
  KEY `variation_set_idx` (`variation_set_id`,`variation_id`)
) ENGINE=MyISAM ;

CREATE TABLE `variation_synonym` (
  `variation_synonym_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `variation_id` int(10) unsigned NOT NULL,
  `subsnp_id` int(15) unsigned DEFAULT NULL,
  `source_id` int(10) unsigned NOT NULL,
  `name` varchar(255) DEFAULT NULL,
  `moltype` varchar(50) DEFAULT NULL,
  PRIMARY KEY (`variation_synonym_id`),
  UNIQUE KEY `name` (`name`,`source_id`),
  KEY `variation_idx` (`variation_id`),
  KEY `subsnp_idx` (`subsnp_id`),
  KEY `source_idx` (`source_id`)
) ENGINE=MyISAM  ;

