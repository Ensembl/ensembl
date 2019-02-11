CREATE TABLE `alt_allele` (
  `alt_allele_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `alt_allele_group_id` int(10) unsigned NOT NULL,
  `gene_id` int(10) unsigned NOT NULL,
  PRIMARY KEY (`alt_allele_id`),
  UNIQUE KEY `gene_idx` (`gene_id`),
  KEY `gene_id` (`gene_id`,`alt_allele_group_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `alt_allele_attrib` (
  `alt_allele_id` int(10) unsigned DEFAULT NULL,
  `attrib` enum('IS_REPRESENTATIVE','IS_MOST_COMMON_ALLELE','IN_CORRECTED_ASSEMBLY','HAS_CODING_POTENTIAL','IN_ARTIFICIALLY_DUPLICATED_ASSEMBLY','IN_SYNTENIC_REGION','HAS_SAME_UNDERLYING_DNA_SEQUENCE','IN_BROKEN_ASSEMBLY_REGION','IS_VALID_ALTERNATE','SAME_AS_REPRESENTATIVE','SAME_AS_ANOTHER_ALLELE','MANUALLY_ASSIGNED','AUTOMATICALLY_ASSIGNED') DEFAULT NULL,
  KEY `aa_idx` (`alt_allele_id`,`attrib`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `alt_allele_group` (
  `alt_allele_group_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  PRIMARY KEY (`alt_allele_group_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `analysis` (
  `analysis_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `created` datetime DEFAULT NULL,
  `logic_name` varchar(40) COLLATE latin1_bin NOT NULL DEFAULT '',
  `db` varchar(120) COLLATE latin1_bin DEFAULT NULL,
  `db_version` varchar(40) COLLATE latin1_bin DEFAULT NULL,
  `db_file` varchar(120) COLLATE latin1_bin DEFAULT NULL,
  `program` varchar(80) COLLATE latin1_bin DEFAULT NULL,
  `program_version` varchar(40) COLLATE latin1_bin DEFAULT NULL,
  `program_file` varchar(80) COLLATE latin1_bin DEFAULT NULL,
  `parameters` text COLLATE latin1_bin,
  `module` varchar(80) COLLATE latin1_bin DEFAULT NULL,
  `module_version` varchar(40) COLLATE latin1_bin DEFAULT NULL,
  `gff_source` varchar(40) COLLATE latin1_bin DEFAULT NULL,
  `gff_feature` varchar(40) COLLATE latin1_bin DEFAULT NULL,
  PRIMARY KEY (`analysis_id`),
  UNIQUE KEY `logic_name_idx` (`logic_name`)
) ENGINE=MyISAM AUTO_INCREMENT=39 DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

CREATE TABLE `analysis_description` (
  `analysis_id` int(10) unsigned NOT NULL DEFAULT '0',
  `description` text COLLATE latin1_bin,
  `display_label` varchar(255) COLLATE latin1_bin DEFAULT NULL,
  `displayable` tinyint(1) NOT NULL DEFAULT '1',
  `web_data` text COLLATE latin1_bin,
  UNIQUE KEY `analysis_idx` (`analysis_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

CREATE TABLE `assembly` (
  `asm_seq_region_id` int(10) unsigned NOT NULL DEFAULT '0',
  `cmp_seq_region_id` int(10) unsigned NOT NULL DEFAULT '0',
  `asm_start` int(10) NOT NULL DEFAULT '0',
  `asm_end` int(10) NOT NULL DEFAULT '0',
  `cmp_start` int(10) NOT NULL DEFAULT '0',
  `cmp_end` int(10) NOT NULL DEFAULT '0',
  `ori` tinyint(4) NOT NULL DEFAULT '0',
  UNIQUE KEY `all_idx` (`asm_seq_region_id`,`cmp_seq_region_id`,`asm_start`,`asm_end`,`cmp_start`,`cmp_end`,`ori`),
  KEY `cmp_seq_region_id` (`cmp_seq_region_id`),
  KEY `asm_seq_region_id` (`asm_seq_region_id`,`asm_start`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

CREATE TABLE `assembly_exception` (
  `assembly_exception_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `seq_region_id` int(11) NOT NULL DEFAULT '0',
  `seq_region_start` int(11) NOT NULL DEFAULT '0',
  `seq_region_end` int(11) NOT NULL DEFAULT '0',
  `exc_type` enum('HAP','PAR','PATCH_NOVEL','PATCH_FIX') COLLATE latin1_bin NOT NULL DEFAULT 'HAP',
  `exc_seq_region_id` int(11) NOT NULL DEFAULT '0',
  `exc_seq_region_start` int(11) NOT NULL DEFAULT '0',
  `exc_seq_region_end` int(11) NOT NULL DEFAULT '0',
  `ori` int(11) NOT NULL DEFAULT '0',
  PRIMARY KEY (`assembly_exception_id`),
  KEY `sr_idx` (`seq_region_id`,`seq_region_start`),
  KEY `ex_idx` (`exc_seq_region_id`,`exc_seq_region_start`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

CREATE TABLE `associated_group` (
  `associated_group_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `description` varchar(128) DEFAULT NULL,
  PRIMARY KEY (`associated_group_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `associated_xref` (
  `associated_xref_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `object_xref_id` int(10) unsigned NOT NULL DEFAULT '0',
  `xref_id` int(10) unsigned NOT NULL DEFAULT '0',
  `source_xref_id` int(10) unsigned DEFAULT NULL,
  `condition_type` varchar(128) DEFAULT NULL,
  `associated_group_id` int(10) unsigned DEFAULT NULL,
  `rank` int(10) unsigned DEFAULT '0',
  PRIMARY KEY (`associated_xref_id`),
  UNIQUE KEY `object_associated_source_type_idx` (`object_xref_id`,`xref_id`,`source_xref_id`,`condition_type`,`associated_group_id`),
  KEY `associated_source_idx` (`source_xref_id`),
  KEY `associated_object_idx` (`object_xref_id`),
  KEY `associated_idx` (`xref_id`),
  KEY `associated_group_idx` (`associated_group_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `attrib_type` (
  `attrib_type_id` smallint(5) unsigned NOT NULL AUTO_INCREMENT,
  `code` varchar(20) COLLATE latin1_bin NOT NULL DEFAULT '',
  `name` varchar(255) COLLATE latin1_bin NOT NULL DEFAULT '',
  `description` text COLLATE latin1_bin,
  PRIMARY KEY (`attrib_type_id`),
  UNIQUE KEY `code_idx` (`code`)
) ENGINE=MyISAM AUTO_INCREMENT=391 DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

CREATE TABLE `biotype` (
  `biotype_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(64) NOT NULL,
  `object_type` enum('gene','transcript') NOT NULL DEFAULT 'gene',
  `db_type` set('cdna','core','coreexpressionatlas','coreexpressionest','coreexpressiongnf','funcgen','otherfeatures','rnaseq','variation','vega','presite','sangervega') NOT NULL DEFAULT 'core',
  `attrib_type_id` int(11) DEFAULT NULL,
  `description` text,
  `biotype_group` enum('coding','pseudogene','snoncoding','lnoncoding','mnoncoding','LRG','undefined','no_group') DEFAULT NULL,
  `so_acc` varchar(64) DEFAULT NULL,
  PRIMARY KEY (`biotype_id`),
  UNIQUE KEY `name_type_idx` (`name`,`object_type`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

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
) ENGINE=MyISAM AUTO_INCREMENT=4 DEFAULT CHARSET=latin1;

CREATE TABLE `data_file` (
  `data_file_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `coord_system_id` int(10) unsigned NOT NULL,
  `analysis_id` smallint(5) unsigned NOT NULL,
  `name` varchar(100) NOT NULL,
  `version_lock` tinyint(1) NOT NULL DEFAULT '0',
  `absolute` tinyint(1) NOT NULL DEFAULT '0',
  `url` text,
  `file_type` enum('BAM','BAMCOV','BIGBED','BIGWIG','VCF') DEFAULT NULL,
  PRIMARY KEY (`data_file_id`),
  UNIQUE KEY `df_unq_idx` (`coord_system_id`,`analysis_id`,`name`,`file_type`),
  KEY `df_name_idx` (`name`),
  KEY `df_analysis_idx` (`analysis_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `density_feature` (
  `density_feature_id` int(11) NOT NULL AUTO_INCREMENT,
  `density_type_id` int(11) NOT NULL DEFAULT '0',
  `seq_region_id` int(11) NOT NULL DEFAULT '0',
  `seq_region_start` int(11) NOT NULL DEFAULT '0',
  `seq_region_end` int(11) NOT NULL DEFAULT '0',
  `density_value` float NOT NULL DEFAULT '0',
  PRIMARY KEY (`density_feature_id`),
  KEY `seq_region_idx` (`density_type_id`,`seq_region_id`,`seq_region_start`),
  KEY `seq_region_id_idx` (`seq_region_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

CREATE TABLE `density_type` (
  `density_type_id` int(11) NOT NULL AUTO_INCREMENT,
  `analysis_id` int(11) NOT NULL DEFAULT '0',
  `block_size` int(11) NOT NULL DEFAULT '0',
  `region_features` int(11) NOT NULL DEFAULT '0',
  `value_type` enum('sum','ratio') COLLATE latin1_bin NOT NULL DEFAULT 'sum',
  PRIMARY KEY (`density_type_id`),
  UNIQUE KEY `analysis_id` (`analysis_id`,`block_size`,`region_features`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

CREATE TABLE `dependent_xref` (
  `object_xref_id` int(10) unsigned NOT NULL,
  `master_xref_id` int(10) unsigned NOT NULL,
  `dependent_xref_id` int(10) unsigned NOT NULL,
  PRIMARY KEY (`object_xref_id`),
  KEY `dependent` (`dependent_xref_id`),
  KEY `master_idx` (`master_xref_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `ditag` (
  `ditag_id` int(10) NOT NULL AUTO_INCREMENT,
  `name` varchar(30) DEFAULT NULL,
  `type` varchar(30) DEFAULT NULL,
  `tag_count` smallint(6) DEFAULT '1',
  `sequence` text,
  PRIMARY KEY (`ditag_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `ditag_feature` (
  `ditag_feature_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `ditag_id` int(10) unsigned NOT NULL DEFAULT '0',
  `ditag_pair_id` int(10) unsigned NOT NULL DEFAULT '0',
  `seq_region_id` int(10) unsigned NOT NULL DEFAULT '0',
  `seq_region_start` int(10) unsigned NOT NULL DEFAULT '0',
  `seq_region_end` int(10) unsigned NOT NULL DEFAULT '0',
  `seq_region_strand` tinyint(1) NOT NULL DEFAULT '0',
  `analysis_id` int(10) unsigned NOT NULL DEFAULT '0',
  `hit_start` int(10) unsigned NOT NULL DEFAULT '0',
  `hit_end` int(10) unsigned NOT NULL DEFAULT '0',
  `hit_strand` tinyint(1) NOT NULL DEFAULT '0',
  `cigar_line` text,
  `ditag_side` char(1) DEFAULT '',
  PRIMARY KEY (`ditag_feature_id`),
  KEY `ditag_id` (`ditag_id`),
  KEY `ditag_pair_id` (`ditag_pair_id`),
  KEY `seq_region_idx` (`seq_region_id`,`seq_region_start`,`seq_region_end`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `dna` (
  `seq_region_id` int(10) unsigned NOT NULL DEFAULT '0',
  `sequence` mediumtext COLLATE latin1_bin NOT NULL,
  PRIMARY KEY (`seq_region_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin MAX_ROWS=750000 AVG_ROW_LENGTH=19000;

CREATE TABLE `dna_align_feature` (
  `dna_align_feature_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `seq_region_id` int(10) unsigned NOT NULL DEFAULT '0',
  `seq_region_start` int(10) unsigned NOT NULL DEFAULT '0',
  `seq_region_end` int(10) unsigned NOT NULL DEFAULT '0',
  `seq_region_strand` tinyint(1) NOT NULL DEFAULT '0',
  `hit_start` int(11) NOT NULL DEFAULT '0',
  `hit_end` int(11) NOT NULL DEFAULT '0',
  `hit_strand` tinyint(1) NOT NULL DEFAULT '0',
  `hit_name` varchar(40) COLLATE latin1_bin NOT NULL DEFAULT '',
  `analysis_id` int(10) unsigned NOT NULL DEFAULT '0',
  `score` double DEFAULT NULL,
  `evalue` double DEFAULT NULL,
  `perc_ident` float DEFAULT NULL,
  `cigar_line` text COLLATE latin1_bin,
  `external_db_id` smallint(5) unsigned DEFAULT NULL,
  `hcoverage` double DEFAULT NULL,
  `align_type` enum('ensembl','cigar','vulgar','mdtag') COLLATE latin1_bin DEFAULT 'ensembl',
  PRIMARY KEY (`dna_align_feature_id`),
  KEY `seq_region_idx` (`seq_region_id`,`analysis_id`,`seq_region_start`,`score`),
  KEY `seq_region_idx_2` (`seq_region_id`,`seq_region_start`),
  KEY `hit_idx` (`hit_name`),
  KEY `analysis_idx` (`analysis_id`),
  KEY `external_db_idx` (`external_db_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin MAX_ROWS=100000000 AVG_ROW_LENGTH=80;

CREATE TABLE `exon` (
  `exon_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `seq_region_id` int(10) unsigned NOT NULL,
  `seq_region_start` int(10) unsigned NOT NULL,
  `seq_region_end` int(10) unsigned NOT NULL,
  `seq_region_strand` tinyint(2) NOT NULL,
  `phase` tinyint(2) NOT NULL,
  `end_phase` tinyint(2) NOT NULL,
  `is_current` tinyint(1) NOT NULL DEFAULT '1',
  `is_constitutive` tinyint(1) NOT NULL DEFAULT '0',
  `stable_id` varchar(128) DEFAULT NULL,
  `version` smallint(5) unsigned DEFAULT NULL,
  `created_date` datetime DEFAULT NULL,
  `modified_date` datetime DEFAULT NULL,
  PRIMARY KEY (`exon_id`),
  KEY `seq_region_idx` (`seq_region_id`,`seq_region_start`),
  KEY `stable_id_idx` (`stable_id`,`version`)
) ENGINE=MyISAM AUTO_INCREMENT=6885 DEFAULT CHARSET=latin1;

CREATE TABLE `exon_transcript` (
  `exon_id` int(10) unsigned NOT NULL DEFAULT '0',
  `transcript_id` int(10) unsigned NOT NULL DEFAULT '0',
  `rank` int(10) NOT NULL DEFAULT '0',
  PRIMARY KEY (`exon_id`,`transcript_id`,`rank`),
  KEY `transcript` (`transcript_id`),
  KEY `exon` (`exon_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

CREATE TABLE `external_db` (
  `external_db_id` int(11) NOT NULL DEFAULT '0',
  `db_name` varchar(27) COLLATE latin1_bin NOT NULL DEFAULT '',
  `db_release` varchar(40) COLLATE latin1_bin NOT NULL DEFAULT '',
  `status` enum('KNOWNXREF','KNOWN','XREF','PRED','ORTH','PSEUDO') COLLATE latin1_bin NOT NULL DEFAULT 'KNOWNXREF',
  `priority` int(11) NOT NULL DEFAULT '0',
  `db_display_name` varchar(255) COLLATE latin1_bin DEFAULT NULL,
  `type` enum('ARRAY','ALT_TRANS','ALT_GENE','MISC','LIT','PRIMARY_DB_SYNONYM','ENSEMBL') COLLATE latin1_bin DEFAULT NULL,
  `secondary_db_name` varchar(255) COLLATE latin1_bin DEFAULT NULL,
  `secondary_db_table` varchar(255) COLLATE latin1_bin DEFAULT NULL,
  `description` text COLLATE latin1_bin,
  PRIMARY KEY (`external_db_id`),
  UNIQUE KEY `db_name_db_release_idx` (`db_name`,`db_release`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

CREATE TABLE `external_synonym` (
  `xref_id` int(10) unsigned NOT NULL DEFAULT '0',
  `synonym` varchar(40) COLLATE latin1_bin NOT NULL DEFAULT '',
  PRIMARY KEY (`xref_id`,`synonym`),
  KEY `name_index` (`synonym`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

CREATE TABLE `gene` (
  `gene_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `biotype` varchar(40) NOT NULL,
  `analysis_id` smallint(5) unsigned NOT NULL,
  `seq_region_id` int(10) unsigned NOT NULL,
  `seq_region_start` int(10) unsigned NOT NULL,
  `seq_region_end` int(10) unsigned NOT NULL,
  `seq_region_strand` tinyint(2) NOT NULL,
  `display_xref_id` int(10) unsigned DEFAULT NULL,
  `source` varchar(40) NOT NULL,
  `description` text,
  `is_current` tinyint(1) NOT NULL DEFAULT '1',
  `canonical_transcript_id` int(10) unsigned NOT NULL,
  `stable_id` varchar(128) DEFAULT NULL,
  `version` smallint(5) unsigned DEFAULT NULL,
  `created_date` datetime DEFAULT NULL,
  `modified_date` datetime DEFAULT NULL,
  PRIMARY KEY (`gene_id`),
  KEY `seq_region_idx` (`seq_region_id`,`seq_region_start`),
  KEY `xref_id_index` (`display_xref_id`),
  KEY `analysis_idx` (`analysis_id`),
  KEY `stable_id_idx` (`stable_id`,`version`),
  KEY `canonical_transcript_id_idx` (`canonical_transcript_id`)
) ENGINE=MyISAM AUTO_INCREMENT=6885 DEFAULT CHARSET=latin1;

CREATE TABLE `gene_archive` (
  `gene_stable_id` varchar(128) COLLATE latin1_bin NOT NULL DEFAULT '',
  `gene_version` smallint(6) NOT NULL DEFAULT '0',
  `transcript_stable_id` varchar(128) COLLATE latin1_bin NOT NULL DEFAULT '',
  `transcript_version` smallint(6) NOT NULL DEFAULT '0',
  `translation_stable_id` varchar(128) COLLATE latin1_bin NOT NULL DEFAULT '',
  `translation_version` smallint(6) NOT NULL DEFAULT '0',
  `peptide_archive_id` int(11) NOT NULL DEFAULT '0',
  `mapping_session_id` int(11) NOT NULL DEFAULT '0',
  KEY `gene_idx` (`gene_stable_id`,`gene_version`),
  KEY `transcript_idx` (`transcript_stable_id`,`transcript_version`),
  KEY `translation_idx` (`translation_stable_id`,`translation_version`),
  KEY `peptide_archive_id_idx` (`peptide_archive_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

CREATE TABLE `gene_attrib` (
  `gene_id` int(10) unsigned NOT NULL DEFAULT '0',
  `attrib_type_id` smallint(5) unsigned NOT NULL DEFAULT '0',
  `value` text COLLATE latin1_bin NOT NULL,
  UNIQUE KEY `gene_attribx` (`gene_id`,`attrib_type_id`,`value`(500)),
  KEY `type_val_idx` (`attrib_type_id`,`value`(40)),
  KEY `val_only_idx` (`value`(40)),
  KEY `gene_idx` (`gene_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

CREATE TABLE `genome_statistics` (
  `genome_statistics_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `statistic` varchar(128) NOT NULL,
  `value` bigint(11) unsigned NOT NULL DEFAULT '0',
  `species_id` int(10) unsigned DEFAULT '1',
  `attrib_type_id` int(10) unsigned DEFAULT NULL,
  `timestamp` datetime DEFAULT NULL,
  PRIMARY KEY (`genome_statistics_id`),
  UNIQUE KEY `stats_uniq` (`statistic`,`attrib_type_id`,`species_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `identity_xref` (
  `object_xref_id` int(10) unsigned NOT NULL DEFAULT '0',
  `xref_identity` int(5) DEFAULT NULL,
  `ensembl_identity` int(5) DEFAULT NULL,
  `xref_start` int(11) DEFAULT NULL,
  `xref_end` int(11) DEFAULT NULL,
  `ensembl_start` int(11) DEFAULT NULL,
  `ensembl_end` int(11) DEFAULT NULL,
  `cigar_line` text COLLATE latin1_bin,
  `score` double DEFAULT NULL,
  `evalue` double DEFAULT NULL,
  PRIMARY KEY (`object_xref_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

CREATE TABLE `interpro` (
  `interpro_ac` varchar(40) COLLATE latin1_bin NOT NULL DEFAULT '',
  `id` varchar(40) COLLATE latin1_bin NOT NULL,
  UNIQUE KEY `accession_idx` (`interpro_ac`,`id`),
  KEY `id_idx` (`id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

CREATE TABLE `intron_supporting_evidence` (
  `intron_supporting_evidence_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `analysis_id` smallint(5) unsigned NOT NULL,
  `seq_region_id` int(10) unsigned NOT NULL,
  `seq_region_start` int(10) unsigned NOT NULL,
  `seq_region_end` int(10) unsigned NOT NULL,
  `seq_region_strand` tinyint(2) NOT NULL,
  `hit_name` varchar(100) NOT NULL,
  `score` decimal(10,3) DEFAULT NULL,
  `score_type` enum('NONE','DEPTH') DEFAULT 'NONE',
  `is_splice_canonical` tinyint(1) NOT NULL DEFAULT '0',
  PRIMARY KEY (`intron_supporting_evidence_id`),
  UNIQUE KEY `analysis_id` (`analysis_id`,`seq_region_id`,`seq_region_start`,`seq_region_end`,`seq_region_strand`,`hit_name`),
  KEY `seq_region_idx` (`seq_region_id`,`seq_region_start`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `karyotype` (
  `karyotype_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `seq_region_id` int(10) unsigned NOT NULL DEFAULT '0',
  `seq_region_start` int(10) NOT NULL DEFAULT '0',
  `seq_region_end` int(10) NOT NULL DEFAULT '0',
  `band` varchar(40) COLLATE latin1_bin DEFAULT NULL,
  `stain` varchar(40) COLLATE latin1_bin DEFAULT NULL,
  PRIMARY KEY (`karyotype_id`),
  KEY `region_band_idx` (`seq_region_id`,`band`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

CREATE TABLE `map` (
  `map_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `map_name` varchar(30) COLLATE latin1_bin NOT NULL DEFAULT '',
  PRIMARY KEY (`map_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

CREATE TABLE `mapping_session` (
  `mapping_session_id` int(11) NOT NULL AUTO_INCREMENT,
  `old_db_name` varchar(80) COLLATE latin1_bin NOT NULL DEFAULT '',
  `new_db_name` varchar(80) COLLATE latin1_bin NOT NULL DEFAULT '',
  `old_release` varchar(5) COLLATE latin1_bin NOT NULL DEFAULT '',
  `new_release` varchar(5) COLLATE latin1_bin NOT NULL DEFAULT '',
  `old_assembly` varchar(20) COLLATE latin1_bin NOT NULL DEFAULT '',
  `new_assembly` varchar(20) COLLATE latin1_bin NOT NULL DEFAULT '',
  `created` datetime DEFAULT NULL,
  PRIMARY KEY (`mapping_session_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

CREATE TABLE `mapping_set` (
  `mapping_set_id` int(10) unsigned NOT NULL,
  `internal_schema_build` varchar(20) NOT NULL,
  `external_schema_build` varchar(20) NOT NULL,
  PRIMARY KEY (`mapping_set_id`),
  UNIQUE KEY `mapping_idx` (`internal_schema_build`,`external_schema_build`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `marker` (
  `marker_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `display_marker_synonym_id` int(10) unsigned DEFAULT NULL,
  `left_primer` varchar(100) COLLATE latin1_bin NOT NULL DEFAULT '',
  `right_primer` varchar(100) COLLATE latin1_bin NOT NULL DEFAULT '',
  `min_primer_dist` int(10) unsigned NOT NULL DEFAULT '0',
  `max_primer_dist` int(10) unsigned NOT NULL DEFAULT '0',
  `priority` int(11) DEFAULT NULL,
  `type` enum('est','microsatellite') COLLATE latin1_bin DEFAULT NULL,
  PRIMARY KEY (`marker_id`),
  KEY `marker_idx` (`marker_id`,`priority`),
  KEY `display_idx` (`display_marker_synonym_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

CREATE TABLE `marker_feature` (
  `marker_feature_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `marker_id` int(10) unsigned NOT NULL DEFAULT '0',
  `seq_region_id` int(10) unsigned NOT NULL DEFAULT '0',
  `seq_region_start` int(10) unsigned NOT NULL DEFAULT '0',
  `seq_region_end` int(10) unsigned NOT NULL DEFAULT '0',
  `analysis_id` int(10) unsigned NOT NULL DEFAULT '0',
  `map_weight` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`marker_feature_id`),
  KEY `seq_region_idx` (`seq_region_id`,`seq_region_start`),
  KEY `analysis_idx` (`analysis_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

CREATE TABLE `marker_map_location` (
  `marker_id` int(10) unsigned NOT NULL DEFAULT '0',
  `map_id` int(10) unsigned NOT NULL DEFAULT '0',
  `chromosome_name` varchar(15) COLLATE latin1_bin NOT NULL DEFAULT '',
  `marker_synonym_id` int(10) unsigned NOT NULL DEFAULT '0',
  `position` varchar(15) COLLATE latin1_bin NOT NULL DEFAULT '',
  `lod_score` double DEFAULT NULL,
  PRIMARY KEY (`marker_id`,`map_id`),
  KEY `map_idx` (`map_id`,`chromosome_name`,`position`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

CREATE TABLE `marker_synonym` (
  `marker_synonym_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `marker_id` int(10) unsigned NOT NULL DEFAULT '0',
  `source` varchar(20) COLLATE latin1_bin DEFAULT NULL,
  `name` varchar(30) COLLATE latin1_bin DEFAULT NULL,
  PRIMARY KEY (`marker_synonym_id`),
  KEY `marker_synonym_idx` (`marker_synonym_id`,`name`),
  KEY `marker_idx` (`marker_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

CREATE TABLE `meta` (
  `meta_id` int(11) NOT NULL AUTO_INCREMENT,
  `species_id` int(10) unsigned DEFAULT '1',
  `meta_key` varchar(40) NOT NULL,
  `meta_value` varchar(255) NOT NULL,
  PRIMARY KEY (`meta_id`),
  UNIQUE KEY `species_key_value_idx` (`species_id`,`meta_key`,`meta_value`),
  KEY `species_value_idx` (`species_id`,`meta_value`)
) ENGINE=MyISAM AUTO_INCREMENT=224 DEFAULT CHARSET=latin1;

CREATE TABLE `meta_coord` (
  `table_name` varchar(40) COLLATE latin1_bin NOT NULL DEFAULT '',
  `coord_system_id` int(11) NOT NULL DEFAULT '0',
  `max_length` int(11) DEFAULT NULL,
  UNIQUE KEY `cs_table_name_idx` (`coord_system_id`,`table_name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

CREATE TABLE `misc_attrib` (
  `misc_feature_id` int(10) unsigned NOT NULL DEFAULT '0',
  `attrib_type_id` smallint(5) unsigned NOT NULL DEFAULT '0',
  `value` text COLLATE latin1_bin NOT NULL,
  UNIQUE KEY `misc_attribx` (`misc_feature_id`,`attrib_type_id`,`value`(500)),
  KEY `type_val_idx` (`attrib_type_id`,`value`(40)),
  KEY `val_only_idx` (`value`(40)),
  KEY `misc_feature_idx` (`misc_feature_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

CREATE TABLE `misc_feature` (
  `misc_feature_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `seq_region_id` int(10) unsigned NOT NULL DEFAULT '0',
  `seq_region_start` int(10) unsigned NOT NULL DEFAULT '0',
  `seq_region_end` int(10) unsigned NOT NULL DEFAULT '0',
  `seq_region_strand` tinyint(4) NOT NULL DEFAULT '0',
  PRIMARY KEY (`misc_feature_id`),
  KEY `seq_region_idx` (`seq_region_id`,`seq_region_start`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

CREATE TABLE `misc_feature_misc_set` (
  `misc_feature_id` int(10) unsigned NOT NULL DEFAULT '0',
  `misc_set_id` smallint(5) unsigned NOT NULL DEFAULT '0',
  PRIMARY KEY (`misc_feature_id`,`misc_set_id`),
  KEY `reverse_idx` (`misc_set_id`,`misc_feature_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

CREATE TABLE `misc_set` (
  `misc_set_id` smallint(5) unsigned NOT NULL AUTO_INCREMENT,
  `code` varchar(25) COLLATE latin1_bin NOT NULL DEFAULT '',
  `name` varchar(255) COLLATE latin1_bin NOT NULL DEFAULT '',
  `description` text COLLATE latin1_bin NOT NULL,
  `max_length` int(10) unsigned NOT NULL DEFAULT '0',
  PRIMARY KEY (`misc_set_id`),
  UNIQUE KEY `code_idx` (`code`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

CREATE TABLE `object_xref` (
  `object_xref_id` int(11) NOT NULL AUTO_INCREMENT,
  `ensembl_id` int(10) unsigned NOT NULL DEFAULT '0',
  `ensembl_object_type` enum('RawContig','Transcript','Gene','Translation','regulatory_factor','regulatory_feature','Marker') COLLATE latin1_bin NOT NULL DEFAULT 'RawContig',
  `xref_id` int(10) unsigned NOT NULL,
  `linkage_annotation` varchar(255) COLLATE latin1_bin DEFAULT NULL,
  `analysis_id` smallint(5) unsigned DEFAULT NULL,
  PRIMARY KEY (`object_xref_id`),
  UNIQUE KEY `xref_idx` (`xref_id`,`ensembl_object_type`,`ensembl_id`),
  KEY `ensembl_idx` (`ensembl_object_type`,`ensembl_id`),
  KEY `analysis_idx` (`analysis_id`)
) ENGINE=MyISAM AUTO_INCREMENT=81424 DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

CREATE TABLE `ontology_xref` (
  `object_xref_id` int(10) unsigned NOT NULL DEFAULT '0',
  `linkage_type` varchar(3) DEFAULT NULL,
  `source_xref_id` int(10) unsigned DEFAULT NULL,
  UNIQUE KEY `object_source_type_idx` (`object_xref_id`,`source_xref_id`,`linkage_type`),
  KEY `object_idx` (`object_xref_id`),
  KEY `source_idx` (`source_xref_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `operon` (
  `operon_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `seq_region_id` int(10) unsigned NOT NULL,
  `seq_region_start` int(10) unsigned NOT NULL,
  `seq_region_end` int(10) unsigned NOT NULL,
  `seq_region_strand` tinyint(2) NOT NULL,
  `display_label` varchar(255) DEFAULT NULL,
  `analysis_id` smallint(5) unsigned NOT NULL,
  `stable_id` varchar(128) DEFAULT NULL,
  `version` smallint(5) unsigned DEFAULT NULL,
  `created_date` datetime DEFAULT NULL,
  `modified_date` datetime DEFAULT NULL,
  PRIMARY KEY (`operon_id`),
  KEY `seq_region_idx` (`seq_region_id`,`seq_region_start`),
  KEY `name_idx` (`display_label`),
  KEY `stable_id_idx` (`stable_id`,`version`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `operon_transcript` (
  `operon_transcript_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `seq_region_id` int(10) unsigned NOT NULL,
  `seq_region_start` int(10) unsigned NOT NULL,
  `seq_region_end` int(10) unsigned NOT NULL,
  `seq_region_strand` tinyint(2) NOT NULL,
  `operon_id` int(10) unsigned NOT NULL,
  `display_label` varchar(255) DEFAULT NULL,
  `analysis_id` smallint(5) unsigned NOT NULL,
  `stable_id` varchar(128) DEFAULT NULL,
  `version` smallint(5) unsigned DEFAULT NULL,
  `created_date` datetime DEFAULT NULL,
  `modified_date` datetime DEFAULT NULL,
  PRIMARY KEY (`operon_transcript_id`),
  KEY `operon_idx` (`operon_id`),
  KEY `seq_region_idx` (`seq_region_id`,`seq_region_start`),
  KEY `stable_id_idx` (`stable_id`,`version`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `operon_transcript_gene` (
  `operon_transcript_id` int(10) unsigned DEFAULT NULL,
  `gene_id` int(10) unsigned DEFAULT NULL,
  KEY `operon_transcript_gene_idx` (`operon_transcript_id`,`gene_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `peptide_archive` (
  `peptide_archive_id` int(11) NOT NULL AUTO_INCREMENT,
  `md5_checksum` varchar(32) COLLATE latin1_bin DEFAULT NULL,
  `peptide_seq` mediumtext COLLATE latin1_bin NOT NULL,
  PRIMARY KEY (`peptide_archive_id`),
  KEY `checksum` (`md5_checksum`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

CREATE TABLE `prediction_exon` (
  `prediction_exon_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `prediction_transcript_id` int(10) unsigned NOT NULL DEFAULT '0',
  `exon_rank` smallint(5) unsigned NOT NULL DEFAULT '0',
  `seq_region_id` int(10) unsigned NOT NULL DEFAULT '0',
  `seq_region_start` int(10) unsigned NOT NULL DEFAULT '0',
  `seq_region_end` int(10) unsigned NOT NULL DEFAULT '0',
  `seq_region_strand` tinyint(4) NOT NULL DEFAULT '0',
  `start_phase` tinyint(4) NOT NULL DEFAULT '0',
  `score` double DEFAULT NULL,
  `p_value` double DEFAULT NULL,
  PRIMARY KEY (`prediction_exon_id`),
  KEY `prediction_transcript_id` (`prediction_transcript_id`),
  KEY `seq_region_id` (`seq_region_id`,`seq_region_start`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

CREATE TABLE `prediction_transcript` (
  `prediction_transcript_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `seq_region_id` int(10) unsigned NOT NULL DEFAULT '0',
  `seq_region_start` int(10) unsigned NOT NULL DEFAULT '0',
  `seq_region_end` int(10) unsigned NOT NULL DEFAULT '0',
  `seq_region_strand` tinyint(4) NOT NULL DEFAULT '0',
  `analysis_id` int(11) DEFAULT NULL,
  `display_label` varchar(255) COLLATE latin1_bin DEFAULT NULL,
  PRIMARY KEY (`prediction_transcript_id`),
  KEY `seq_region_idx` (`seq_region_id`,`seq_region_start`),
  KEY `analysis_idx` (`analysis_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

CREATE TABLE `protein_align_feature` (
  `protein_align_feature_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `seq_region_id` int(10) unsigned NOT NULL DEFAULT '0',
  `seq_region_start` int(10) unsigned NOT NULL DEFAULT '0',
  `seq_region_end` int(10) unsigned NOT NULL DEFAULT '0',
  `seq_region_strand` tinyint(1) NOT NULL DEFAULT '1',
  `hit_start` int(10) NOT NULL DEFAULT '0',
  `hit_end` int(10) NOT NULL DEFAULT '0',
  `hit_name` varchar(40) COLLATE latin1_bin NOT NULL DEFAULT '',
  `analysis_id` int(10) unsigned NOT NULL DEFAULT '0',
  `score` double DEFAULT NULL,
  `evalue` double DEFAULT NULL,
  `perc_ident` float DEFAULT NULL,
  `cigar_line` text COLLATE latin1_bin,
  `external_db_id` smallint(5) unsigned DEFAULT NULL,
  `hcoverage` double DEFAULT NULL,
  `align_type` enum('ensembl','cigar','vulgar','mdtag') COLLATE latin1_bin DEFAULT 'ensembl',
  PRIMARY KEY (`protein_align_feature_id`),
  KEY `seq_region_idx` (`seq_region_id`,`analysis_id`,`seq_region_start`,`score`),
  KEY `seq_region_idx_2` (`seq_region_id`,`seq_region_start`),
  KEY `hit_idx` (`hit_name`),
  KEY `analysis_idx` (`analysis_id`),
  KEY `external_db_idx` (`external_db_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin MAX_ROWS=100000000 AVG_ROW_LENGTH=80;

CREATE TABLE `protein_feature` (
  `protein_feature_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `translation_id` int(11) NOT NULL DEFAULT '0',
  `seq_start` int(10) NOT NULL DEFAULT '0',
  `seq_end` int(10) NOT NULL DEFAULT '0',
  `hit_start` int(10) NOT NULL DEFAULT '0',
  `hit_end` int(10) NOT NULL DEFAULT '0',
  `hit_name` varchar(40) COLLATE latin1_bin NOT NULL,
  `analysis_id` int(10) unsigned NOT NULL DEFAULT '0',
  `score` double NOT NULL DEFAULT '0',
  `evalue` double DEFAULT NULL,
  `perc_ident` float DEFAULT NULL,
  `external_data` text COLLATE latin1_bin,
  `hit_description` text COLLATE latin1_bin,
  `cigar_line` text COLLATE latin1_bin,
  `align_type` enum('ensembl','cigar','cigarplus','vulgar','mdtag') COLLATE latin1_bin DEFAULT NULL,
  PRIMARY KEY (`protein_feature_id`),
  UNIQUE KEY `aln_idx` (`translation_id`,`hit_name`,`seq_start`,`seq_end`,`hit_start`,`hit_end`,`analysis_id`),
  KEY `translation_idx` (`translation_id`),
  KEY `hitname_idx` (`hit_name`),
  KEY `analysis_idx` (`analysis_id`)
) ENGINE=MyISAM AUTO_INCREMENT=24117 DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

CREATE TABLE `repeat_consensus` (
  `repeat_consensus_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `repeat_name` varchar(255) COLLATE latin1_bin NOT NULL DEFAULT '',
  `repeat_class` varchar(100) COLLATE latin1_bin NOT NULL DEFAULT '',
  `repeat_type` varchar(40) COLLATE latin1_bin NOT NULL DEFAULT '',
  `repeat_consensus` text COLLATE latin1_bin,
  PRIMARY KEY (`repeat_consensus_id`),
  KEY `name` (`repeat_name`),
  KEY `class` (`repeat_class`),
  KEY `consensus` (`repeat_consensus`(10)),
  KEY `type` (`repeat_type`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

CREATE TABLE `repeat_feature` (
  `repeat_feature_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `seq_region_id` int(10) unsigned NOT NULL DEFAULT '0',
  `seq_region_start` int(10) unsigned NOT NULL DEFAULT '0',
  `seq_region_end` int(10) unsigned NOT NULL DEFAULT '0',
  `seq_region_strand` tinyint(1) NOT NULL DEFAULT '1',
  `repeat_start` int(10) NOT NULL DEFAULT '0',
  `repeat_end` int(10) NOT NULL DEFAULT '0',
  `repeat_consensus_id` int(10) unsigned NOT NULL DEFAULT '0',
  `analysis_id` int(10) unsigned NOT NULL DEFAULT '0',
  `score` double DEFAULT NULL,
  PRIMARY KEY (`repeat_feature_id`),
  KEY `seq_region_idx` (`seq_region_id`,`seq_region_start`),
  KEY `repeat_idx` (`repeat_consensus_id`),
  KEY `analysis_idx` (`analysis_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin MAX_ROWS=100000000 AVG_ROW_LENGTH=80;

CREATE TABLE `seq_region` (
  `seq_region_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(255) COLLATE latin1_bin NOT NULL,
  `coord_system_id` int(10) NOT NULL DEFAULT '0',
  `length` int(10) NOT NULL DEFAULT '0',
  PRIMARY KEY (`seq_region_id`),
  UNIQUE KEY `name_cs_idx` (`name`,`coord_system_id`),
  KEY `cs_idx` (`coord_system_id`)
) ENGINE=MyISAM AUTO_INCREMENT=14 DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

CREATE TABLE `seq_region_attrib` (
  `seq_region_id` int(10) unsigned NOT NULL DEFAULT '0',
  `attrib_type_id` smallint(5) unsigned NOT NULL DEFAULT '0',
  `value` text COLLATE latin1_bin NOT NULL,
  UNIQUE KEY `region_attribx` (`seq_region_id`,`attrib_type_id`,`value`(500)),
  KEY `type_val_idx` (`attrib_type_id`,`value`(40)),
  KEY `val_only_idx` (`value`(40)),
  KEY `seq_region_idx` (`seq_region_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

CREATE TABLE `seq_region_mapping` (
  `external_seq_region_id` int(10) unsigned NOT NULL,
  `internal_seq_region_id` int(10) unsigned NOT NULL,
  `mapping_set_id` int(10) unsigned NOT NULL,
  KEY `mapping_set_idx` (`mapping_set_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `seq_region_synonym` (
  `seq_region_synonym_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `seq_region_id` int(10) unsigned NOT NULL,
  `synonym` varchar(250) NOT NULL,
  `external_db_id` smallint(5) unsigned DEFAULT NULL,
  PRIMARY KEY (`seq_region_synonym_id`),
  UNIQUE KEY `syn_idx` (`synonym`,`seq_region_id`),
  KEY `seq_region_idx` (`seq_region_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `simple_feature` (
  `simple_feature_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `seq_region_id` int(10) unsigned NOT NULL DEFAULT '0',
  `seq_region_start` int(10) unsigned NOT NULL DEFAULT '0',
  `seq_region_end` int(10) unsigned NOT NULL DEFAULT '0',
  `seq_region_strand` tinyint(1) NOT NULL DEFAULT '0',
  `display_label` varchar(40) COLLATE latin1_bin NOT NULL DEFAULT '',
  `analysis_id` int(10) unsigned NOT NULL DEFAULT '0',
  `score` double DEFAULT NULL,
  PRIMARY KEY (`simple_feature_id`),
  KEY `seq_region_idx` (`seq_region_id`,`seq_region_start`),
  KEY `analysis_idx` (`analysis_id`),
  KEY `hit_idx` (`display_label`)
) ENGINE=MyISAM AUTO_INCREMENT=37 DEFAULT CHARSET=latin1 COLLATE=latin1_bin MAX_ROWS=100000000 AVG_ROW_LENGTH=80;

CREATE TABLE `stable_id_event` (
  `old_stable_id` varchar(128) COLLATE latin1_bin DEFAULT NULL,
  `old_version` smallint(6) DEFAULT NULL,
  `new_stable_id` varchar(128) COLLATE latin1_bin DEFAULT NULL,
  `new_version` smallint(6) DEFAULT NULL,
  `mapping_session_id` int(10) NOT NULL DEFAULT '0',
  `type` enum('gene','transcript','translation') COLLATE latin1_bin NOT NULL DEFAULT 'gene',
  `score` float NOT NULL DEFAULT '0',
  UNIQUE KEY `uni_idx` (`mapping_session_id`,`old_stable_id`,`old_version`,`new_stable_id`,`new_version`,`type`),
  KEY `new_idx` (`new_stable_id`),
  KEY `old_idx` (`old_stable_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

CREATE TABLE `supporting_feature` (
  `exon_id` int(11) NOT NULL DEFAULT '0',
  `feature_type` enum('dna_align_feature','protein_align_feature') COLLATE latin1_bin DEFAULT NULL,
  `feature_id` int(11) NOT NULL DEFAULT '0',
  UNIQUE KEY `all_idx` (`exon_id`,`feature_type`,`feature_id`),
  KEY `feature_idx` (`feature_type`,`feature_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin MAX_ROWS=100000000 AVG_ROW_LENGTH=80;

CREATE TABLE `transcript` (
  `transcript_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `gene_id` int(10) unsigned DEFAULT NULL,
  `analysis_id` smallint(5) unsigned NOT NULL,
  `seq_region_id` int(10) unsigned NOT NULL,
  `seq_region_start` int(10) unsigned NOT NULL,
  `seq_region_end` int(10) unsigned NOT NULL,
  `seq_region_strand` tinyint(2) NOT NULL,
  `display_xref_id` int(10) unsigned DEFAULT NULL,
  `source` varchar(40) NOT NULL DEFAULT 'ensembl',
  `biotype` varchar(40) NOT NULL,
  `description` text,
  `is_current` tinyint(1) NOT NULL DEFAULT '1',
  `canonical_translation_id` int(10) unsigned DEFAULT NULL,
  `stable_id` varchar(128) DEFAULT NULL,
  `version` smallint(5) unsigned DEFAULT NULL,
  `created_date` datetime DEFAULT NULL,
  `modified_date` datetime DEFAULT NULL,
  PRIMARY KEY (`transcript_id`),
  UNIQUE KEY `canonical_translation_idx` (`canonical_translation_id`),
  KEY `seq_region_idx` (`seq_region_id`,`seq_region_start`),
  KEY `gene_index` (`gene_id`),
  KEY `xref_id_index` (`display_xref_id`),
  KEY `analysis_idx` (`analysis_id`),
  KEY `stable_id_idx` (`stable_id`,`version`)
) ENGINE=MyISAM AUTO_INCREMENT=6885 DEFAULT CHARSET=latin1;

CREATE TABLE `transcript_attrib` (
  `transcript_id` int(10) unsigned NOT NULL DEFAULT '0',
  `attrib_type_id` smallint(5) unsigned NOT NULL DEFAULT '0',
  `value` text COLLATE latin1_bin NOT NULL,
  UNIQUE KEY `transcript_attribx` (`transcript_id`,`attrib_type_id`,`value`(500)),
  KEY `type_val_idx` (`attrib_type_id`,`value`(40)),
  KEY `val_only_idx` (`value`(40)),
  KEY `transcript_idx` (`transcript_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

CREATE TABLE `transcript_intron_supporting_evidence` (
  `transcript_id` int(10) unsigned NOT NULL,
  `intron_supporting_evidence_id` int(10) unsigned NOT NULL,
  `previous_exon_id` int(10) unsigned NOT NULL,
  `next_exon_id` int(10) unsigned NOT NULL,
  PRIMARY KEY (`intron_supporting_evidence_id`,`transcript_id`),
  KEY `transcript_idx` (`transcript_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `transcript_supporting_feature` (
  `transcript_id` int(11) NOT NULL DEFAULT '0',
  `feature_type` enum('dna_align_feature','protein_align_feature') COLLATE latin1_bin DEFAULT NULL,
  `feature_id` int(11) NOT NULL DEFAULT '0',
  UNIQUE KEY `all_idx` (`transcript_id`,`feature_type`,`feature_id`),
  KEY `feature_idx` (`feature_type`,`feature_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin MAX_ROWS=100000000 AVG_ROW_LENGTH=80;

CREATE TABLE `translation` (
  `translation_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `transcript_id` int(10) unsigned NOT NULL,
  `seq_start` int(10) NOT NULL,
  `start_exon_id` int(10) unsigned NOT NULL,
  `seq_end` int(10) NOT NULL,
  `end_exon_id` int(10) unsigned NOT NULL,
  `stable_id` varchar(128) DEFAULT NULL,
  `version` smallint(5) unsigned DEFAULT NULL,
  `created_date` datetime DEFAULT NULL,
  `modified_date` datetime DEFAULT NULL,
  PRIMARY KEY (`translation_id`),
  KEY `transcript_idx` (`transcript_id`),
  KEY `stable_id_idx` (`stable_id`,`version`)
) ENGINE=MyISAM AUTO_INCREMENT=6690 DEFAULT CHARSET=latin1;

CREATE TABLE `translation_attrib` (
  `translation_id` int(10) unsigned NOT NULL DEFAULT '0',
  `attrib_type_id` smallint(5) unsigned NOT NULL DEFAULT '0',
  `value` text COLLATE latin1_bin NOT NULL,
  UNIQUE KEY `translation_attribx` (`translation_id`,`attrib_type_id`,`value`(500)),
  KEY `type_val_idx` (`attrib_type_id`,`value`(40)),
  KEY `val_only_idx` (`value`(40)),
  KEY `translation_idx` (`translation_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

CREATE TABLE `unmapped_object` (
  `unmapped_object_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `type` enum('xref','cDNA','Marker') COLLATE latin1_bin NOT NULL,
  `analysis_id` int(10) unsigned NOT NULL,
  `external_db_id` int(11) DEFAULT NULL,
  `identifier` varchar(255) COLLATE latin1_bin NOT NULL,
  `unmapped_reason_id` int(10) unsigned NOT NULL,
  `query_score` double DEFAULT NULL,
  `target_score` double DEFAULT NULL,
  `ensembl_id` int(10) unsigned DEFAULT '0',
  `ensembl_object_type` enum('RawContig','Transcript','Gene','Translation') COLLATE latin1_bin DEFAULT 'RawContig',
  `parent` varchar(255) COLLATE latin1_bin DEFAULT NULL,
  PRIMARY KEY (`unmapped_object_id`),
  UNIQUE KEY `unique_unmapped_obj_idx` (`ensembl_id`,`ensembl_object_type`,`identifier`,`unmapped_reason_id`,`parent`,`external_db_id`),
  KEY `id_idx` (`identifier`(50)),
  KEY `ext_db_identifier_idx` (`external_db_id`,`identifier`),
  KEY `anal_exdb_idx` (`analysis_id`,`external_db_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

CREATE TABLE `unmapped_reason` (
  `unmapped_reason_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `summary_description` varchar(255) COLLATE latin1_bin DEFAULT NULL,
  `full_description` varchar(255) COLLATE latin1_bin DEFAULT NULL,
  PRIMARY KEY (`unmapped_reason_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

CREATE TABLE `xref` (
  `xref_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `external_db_id` int(11) NOT NULL,
  `dbprimary_acc` varchar(512) COLLATE latin1_bin NOT NULL,
  `display_label` varchar(512) COLLATE latin1_bin NOT NULL,
  `version` varchar(10) COLLATE latin1_bin DEFAULT NULL,
  `description` text COLLATE latin1_bin,
  `info_type` enum('NONE','PROJECTION','MISC','DEPENDENT','DIRECT','SEQUENCE_MATCH','INFERRED_PAIR','PROBE','UNMAPPED','CHECKSUM') COLLATE latin1_bin NOT NULL DEFAULT 'NONE',
  `info_text` varchar(255) COLLATE latin1_bin NOT NULL DEFAULT '',
  PRIMARY KEY (`xref_id`),
  UNIQUE KEY `id_index` (`dbprimary_acc`,`external_db_id`,`info_type`,`info_text`,`version`),
  KEY `display_index` (`display_label`),
  KEY `info_type_idx` (`info_type`)
) ENGINE=MyISAM AUTO_INCREMENT=91703 DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

