CREATE TABLE `alt_allele` (
  `alt_allele_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `alt_allele_group_id` integer  NOT NULL
,  `gene_id` integer  NOT NULL
,  UNIQUE (`gene_id`)
);
CREATE TABLE sqlite_sequence(name,seq);
CREATE TABLE `alt_allele_attrib` (
  `alt_allele_id` integer  DEFAULT NULL
,  `attrib` text  DEFAULT NULL
);
CREATE TABLE `alt_allele_group` (
  `alt_allele_group_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
);
CREATE TABLE `analysis` (
  `analysis_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `created` datetime DEFAULT NULL
,  `logic_name` varchar(128) NOT NULL
,  `db` varchar(120) DEFAULT NULL
,  `db_version` varchar(40) DEFAULT NULL
,  `db_file` varchar(120) DEFAULT NULL
,  `program` varchar(80) DEFAULT NULL
,  `program_version` varchar(40) DEFAULT NULL
,  `program_file` varchar(80) DEFAULT NULL
,  `parameters` text
,  `module` varchar(80) DEFAULT NULL
,  `module_version` varchar(40) DEFAULT NULL
,  `gff_source` varchar(40) DEFAULT NULL
,  `gff_feature` varchar(40) DEFAULT NULL
,  UNIQUE (`logic_name`)
);
CREATE TABLE `analysis_description` (
  `analysis_id` integer  NOT NULL
,  `description` text
,  `display_label` varchar(255) NOT NULL
,  `displayable` integer NOT NULL DEFAULT '1'
,  `web_data` text
,  UNIQUE (`analysis_id`)
);
CREATE TABLE `assembly` (
  `asm_seq_region_id` integer  NOT NULL
,  `cmp_seq_region_id` integer  NOT NULL
,  `asm_start` integer NOT NULL
,  `asm_end` integer NOT NULL
,  `cmp_start` integer NOT NULL
,  `cmp_end` integer NOT NULL
,  `ori` integer NOT NULL
,  UNIQUE (`asm_seq_region_id`,`cmp_seq_region_id`,`asm_start`,`asm_end`,`cmp_start`,`cmp_end`,`ori`)
);
CREATE TABLE `assembly_exception` (
  `assembly_exception_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `seq_region_id` integer  NOT NULL
,  `seq_region_start` integer  NOT NULL
,  `seq_region_end` integer  NOT NULL
,  `exc_type` text  NOT NULL
,  `exc_seq_region_id` integer  NOT NULL
,  `exc_seq_region_start` integer  NOT NULL
,  `exc_seq_region_end` integer  NOT NULL
,  `ori` integer NOT NULL
);
CREATE TABLE `associated_group` (
  `associated_group_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `description` varchar(128) DEFAULT NULL
);
CREATE TABLE `associated_xref` (
  `associated_xref_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `object_xref_id` integer  NOT NULL DEFAULT '0'
,  `xref_id` integer  NOT NULL DEFAULT '0'
,  `source_xref_id` integer  DEFAULT NULL
,  `condition_type` varchar(128) DEFAULT NULL
,  `associated_group_id` integer  DEFAULT NULL
,  `rank` integer  DEFAULT '0'
,  UNIQUE (`object_xref_id`,`xref_id`,`source_xref_id`,`condition_type`,`associated_group_id`)
);
CREATE TABLE `attrib_type` (
  `attrib_type_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `code` varchar(20) NOT NULL DEFAULT ''
,  `name` varchar(255) NOT NULL DEFAULT ''
,  `description` text
,  UNIQUE (`code`)
);
CREATE TABLE `biotype` (
  `biotype_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `name` varchar(64) NOT NULL
,  `object_type` text  NOT NULL DEFAULT 'gene'
,  `db_type` text  NOT NULL DEFAULT 'core'
,  `attrib_type_id` integer DEFAULT NULL
,  `description` text
,  `biotype_group` text  DEFAULT NULL
,  `so_acc` varchar(64) DEFAULT NULL
,  `so_term` varchar(1023) DEFAULT NULL
,  UNIQUE (`name`,`object_type`)
);
CREATE TABLE `coord_system` (
  `coord_system_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `species_id` integer  NOT NULL DEFAULT '1'
,  `name` varchar(40) NOT NULL
,  `version` varchar(255) DEFAULT NULL
,  `rank` integer NOT NULL
,  `attrib` text  DEFAULT NULL
,  UNIQUE (`rank`,`species_id`)
,  UNIQUE (`name`,`version`,`species_id`)
);
CREATE TABLE `data_file` (
  `data_file_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `coord_system_id` integer  NOT NULL
,  `analysis_id` integer  NOT NULL
,  `name` varchar(100) NOT NULL
,  `version_lock` integer NOT NULL DEFAULT '0'
,  `absolute` integer NOT NULL DEFAULT '0'
,  `url` text
,  `file_type` text  DEFAULT NULL
,  UNIQUE (`coord_system_id`,`analysis_id`,`name`,`file_type`)
);
CREATE TABLE `density_feature` (
  `density_feature_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `density_type_id` integer  NOT NULL
,  `seq_region_id` integer  NOT NULL
,  `seq_region_start` integer  NOT NULL
,  `seq_region_end` integer  NOT NULL
,  `density_value` float NOT NULL
);
CREATE TABLE `density_type` (
  `density_type_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `analysis_id` integer  NOT NULL
,  `block_size` integer NOT NULL
,  `region_features` integer NOT NULL
,  `value_type` text  NOT NULL
,  UNIQUE (`analysis_id`,`block_size`,`region_features`)
);
CREATE TABLE `dependent_xref` (
  `object_xref_id` integer  NOT NULL
,  `master_xref_id` integer  NOT NULL
,  `dependent_xref_id` integer  NOT NULL
,  PRIMARY KEY (`object_xref_id`)
);
CREATE TABLE `ditag` (
  `ditag_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `name` varchar(30) NOT NULL
,  `type` varchar(30) NOT NULL
,  `tag_count` integer  NOT NULL DEFAULT '1'
,  `sequence` tinytext NOT NULL
);
CREATE TABLE `ditag_feature` (
  `ditag_feature_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `ditag_id` integer  NOT NULL DEFAULT '0'
,  `ditag_pair_id` integer  NOT NULL DEFAULT '0'
,  `seq_region_id` integer  NOT NULL DEFAULT '0'
,  `seq_region_start` integer  NOT NULL DEFAULT '0'
,  `seq_region_end` integer  NOT NULL DEFAULT '0'
,  `seq_region_strand` integer NOT NULL DEFAULT '0'
,  `analysis_id` integer  NOT NULL DEFAULT '0'
,  `hit_start` integer  NOT NULL DEFAULT '0'
,  `hit_end` integer  NOT NULL DEFAULT '0'
,  `hit_strand` integer NOT NULL DEFAULT '0'
,  `cigar_line` tinytext NOT NULL
,  `ditag_side` text  NOT NULL
);
CREATE TABLE `dna` (
  `seq_region_id` integer  NOT NULL
,  `sequence` longtext NOT NULL
,  PRIMARY KEY (`seq_region_id`)
);
CREATE TABLE `dna_align_feature` (
  `dna_align_feature_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `seq_region_id` integer  NOT NULL
,  `seq_region_start` integer  NOT NULL
,  `seq_region_end` integer  NOT NULL
,  `seq_region_strand` integer NOT NULL
,  `hit_start` integer NOT NULL
,  `hit_end` integer NOT NULL
,  `hit_strand` integer NOT NULL
,  `hit_name` varchar(40) NOT NULL
,  `analysis_id` integer  NOT NULL
,  `score` double DEFAULT NULL
,  `evalue` double DEFAULT NULL
,  `perc_ident` float DEFAULT NULL
,  `cigar_line` text
,  `external_db_id` integer  DEFAULT NULL
,  `hcoverage` double DEFAULT NULL
,  `align_type` text  DEFAULT 'ensembl'
);
CREATE TABLE `dna_align_feature_attrib` (
  `dna_align_feature_id` integer  NOT NULL
,  `attrib_type_id` integer  NOT NULL
,  `value` text NOT NULL
,  UNIQUE (`dna_align_feature_id`,`attrib_type_id`,`value`)
);
CREATE TABLE `exon` (
  `exon_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `seq_region_id` integer  NOT NULL
,  `seq_region_start` integer  NOT NULL
,  `seq_region_end` integer  NOT NULL
,  `seq_region_strand` integer NOT NULL
,  `phase` integer NOT NULL
,  `end_phase` integer NOT NULL
,  `is_current` integer NOT NULL DEFAULT '1'
,  `is_constitutive` integer NOT NULL DEFAULT '0'
,  `stable_id` varchar(128) DEFAULT NULL
,  `version` integer  DEFAULT NULL
,  `created_date` datetime DEFAULT NULL
,  `modified_date` datetime DEFAULT NULL
);
CREATE TABLE `exon_transcript` (
  `exon_id` integer  NOT NULL
,  `transcript_id` integer  NOT NULL
,  `rank` integer NOT NULL
,  PRIMARY KEY (`exon_id`,`transcript_id`,`rank`)
);
CREATE TABLE `external_db` (
  `external_db_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `db_name` varchar(100) NOT NULL
,  `db_release` varchar(255) DEFAULT NULL
,  `status` text  NOT NULL
,  `priority` integer NOT NULL
,  `db_display_name` varchar(255) DEFAULT NULL
,  `type` text  NOT NULL
,  `secondary_db_name` varchar(255) DEFAULT NULL
,  `secondary_db_table` varchar(255) DEFAULT NULL
,  `description` text
,  UNIQUE (`db_name`,`db_release`)
);
CREATE TABLE `external_synonym` (
  `xref_id` integer  NOT NULL
,  `synonym` varchar(100) NOT NULL
,  PRIMARY KEY (`xref_id`,`synonym`)
);
CREATE TABLE `gene` (
  `gene_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `biotype` varchar(40) NOT NULL
,  `analysis_id` integer  NOT NULL
,  `seq_region_id` integer  NOT NULL
,  `seq_region_start` integer  NOT NULL
,  `seq_region_end` integer  NOT NULL
,  `seq_region_strand` integer NOT NULL
,  `display_xref_id` integer  DEFAULT NULL
,  `source` varchar(40) NOT NULL
,  `description` text
,  `is_current` integer NOT NULL DEFAULT '1'
,  `canonical_transcript_id` integer  NOT NULL
,  `stable_id` varchar(128) DEFAULT NULL
,  `version` integer  DEFAULT NULL
,  `created_date` datetime DEFAULT NULL
,  `modified_date` datetime DEFAULT NULL
);
CREATE TABLE `gene_archive` (
  `gene_stable_id` varchar(128) NOT NULL
,  `gene_version` integer NOT NULL DEFAULT '1'
,  `transcript_stable_id` varchar(128) NOT NULL
,  `transcript_version` integer NOT NULL DEFAULT '1'
,  `translation_stable_id` varchar(128) DEFAULT NULL
,  `translation_version` integer NOT NULL DEFAULT '1'
,  `peptide_archive_id` integer  DEFAULT NULL
,  `mapping_session_id` integer  NOT NULL
);
CREATE TABLE `gene_attrib` (
  `gene_id` integer  NOT NULL DEFAULT '0'
,  `attrib_type_id` integer  NOT NULL DEFAULT '0'
,  `value` text NOT NULL
,  UNIQUE (`gene_id`,`attrib_type_id`,`value`)
);
CREATE TABLE `genome_statistics` (
  `genome_statistics_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `statistic` varchar(128) NOT NULL
,  `value` integer  NOT NULL DEFAULT '0'
,  `species_id` integer  DEFAULT '1'
,  `attrib_type_id` integer  DEFAULT NULL
,  `timestamp` datetime DEFAULT NULL
,  UNIQUE (`statistic`,`attrib_type_id`,`species_id`)
);
CREATE TABLE `identity_xref` (
  `object_xref_id` integer  NOT NULL
,  `xref_identity` integer DEFAULT NULL
,  `ensembl_identity` integer DEFAULT NULL
,  `xref_start` integer DEFAULT NULL
,  `xref_end` integer DEFAULT NULL
,  `ensembl_start` integer DEFAULT NULL
,  `ensembl_end` integer DEFAULT NULL
,  `cigar_line` text
,  `score` double DEFAULT NULL
,  `evalue` double DEFAULT NULL
,  PRIMARY KEY (`object_xref_id`)
);
CREATE TABLE `interpro` (
  `interpro_ac` varchar(40) NOT NULL
,  `id` varchar(40) NOT NULL
,  UNIQUE (`interpro_ac`,`id`)
);
CREATE TABLE `intron_supporting_evidence` (
  `intron_supporting_evidence_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `analysis_id` integer  NOT NULL
,  `seq_region_id` integer  NOT NULL
,  `seq_region_start` integer  NOT NULL
,  `seq_region_end` integer  NOT NULL
,  `seq_region_strand` integer NOT NULL
,  `hit_name` varchar(100) NOT NULL
,  `score` decimal(10,3) DEFAULT NULL
,  `score_type` text  DEFAULT 'NONE'
,  `is_splice_canonical` integer NOT NULL DEFAULT '0'
,  UNIQUE (`analysis_id`,`seq_region_id`,`seq_region_start`,`seq_region_end`,`seq_region_strand`,`hit_name`)
);
CREATE TABLE `karyotype` (
  `karyotype_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `seq_region_id` integer  NOT NULL
,  `seq_region_start` integer  NOT NULL
,  `seq_region_end` integer  NOT NULL
,  `band` varchar(40) DEFAULT NULL
,  `stain` varchar(40) DEFAULT NULL
);
CREATE TABLE `map` (
  `map_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `map_name` varchar(30) NOT NULL
);
CREATE TABLE `mapping_session` (
  `mapping_session_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `old_db_name` varchar(80) NOT NULL DEFAULT ''
,  `new_db_name` varchar(80) NOT NULL DEFAULT ''
,  `old_release` varchar(5) NOT NULL DEFAULT ''
,  `new_release` varchar(5) NOT NULL DEFAULT ''
,  `old_assembly` varchar(80) NOT NULL DEFAULT ''
,  `new_assembly` varchar(80) NOT NULL DEFAULT ''
,  `created` datetime NOT NULL
);
CREATE TABLE `mapping_set` (
  `mapping_set_id` integer  NOT NULL
,  `internal_schema_build` varchar(20) NOT NULL
,  `external_schema_build` varchar(20) NOT NULL
,  PRIMARY KEY (`mapping_set_id`)
,  UNIQUE (`internal_schema_build`,`external_schema_build`)
);
CREATE TABLE `marker` (
  `marker_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `display_marker_synonym_id` integer  DEFAULT NULL
,  `left_primer` varchar(100) NOT NULL
,  `right_primer` varchar(100) NOT NULL
,  `min_primer_dist` integer  NOT NULL
,  `max_primer_dist` integer  NOT NULL
,  `priority` integer DEFAULT NULL
,  `type` text  DEFAULT NULL
);
CREATE TABLE `marker_feature` (
  `marker_feature_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `marker_id` integer  NOT NULL
,  `seq_region_id` integer  NOT NULL
,  `seq_region_start` integer  NOT NULL
,  `seq_region_end` integer  NOT NULL
,  `analysis_id` integer  NOT NULL
,  `map_weight` integer  DEFAULT NULL
);
CREATE TABLE `marker_map_location` (
  `marker_id` integer  NOT NULL
,  `map_id` integer  NOT NULL
,  `chromosome_name` varchar(15) NOT NULL
,  `marker_synonym_id` integer  NOT NULL
,  `position` varchar(15) NOT NULL
,  `lod_score` double DEFAULT NULL
,  PRIMARY KEY (`marker_id`,`map_id`)
);
CREATE TABLE `marker_synonym` (
  `marker_synonym_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `marker_id` integer  NOT NULL
,  `source` varchar(20) DEFAULT NULL
,  `name` varchar(50) DEFAULT NULL
);
CREATE TABLE `meta` (
  `meta_id` integer NOT NULL PRIMARY KEY AUTOINCREMENT
,  `species_id` integer  DEFAULT '1'
,  `meta_key` varchar(40) NOT NULL
,  `meta_value` varchar(255) NOT NULL
,  UNIQUE (`species_id`,`meta_key`,`meta_value`)
);
CREATE TABLE `meta_coord` (
  `table_name` varchar(40) NOT NULL
,  `coord_system_id` integer  NOT NULL
,  `max_length` integer DEFAULT NULL
,  UNIQUE (`coord_system_id`,`table_name`)
);
CREATE TABLE `misc_attrib` (
  `misc_feature_id` integer  NOT NULL DEFAULT '0'
,  `attrib_type_id` integer  NOT NULL DEFAULT '0'
,  `value` text NOT NULL
,  UNIQUE (`misc_feature_id`,`attrib_type_id`,`value`)
);
CREATE TABLE `misc_feature` (
  `misc_feature_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `seq_region_id` integer  NOT NULL DEFAULT '0'
,  `seq_region_start` integer  NOT NULL DEFAULT '0'
,  `seq_region_end` integer  NOT NULL DEFAULT '0'
,  `seq_region_strand` integer NOT NULL DEFAULT '0'
);
CREATE TABLE `misc_feature_misc_set` (
  `misc_feature_id` integer  NOT NULL DEFAULT '0'
,  `misc_set_id` integer  NOT NULL DEFAULT '0'
,  PRIMARY KEY (`misc_feature_id`,`misc_set_id`)
);
CREATE TABLE `misc_set` (
  `misc_set_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `code` varchar(25) NOT NULL DEFAULT ''
,  `name` varchar(255) NOT NULL DEFAULT ''
,  `description` text NOT NULL
,  `max_length` integer  NOT NULL
,  UNIQUE (`code`)
);
CREATE TABLE `object_xref` (
  `object_xref_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `ensembl_id` integer  NOT NULL
,  `ensembl_object_type` text  NOT NULL
,  `xref_id` integer  NOT NULL
,  `linkage_annotation` varchar(255) DEFAULT NULL
,  `analysis_id` integer  DEFAULT NULL
,  UNIQUE (`xref_id`,`ensembl_object_type`,`ensembl_id`,`analysis_id`)
);
CREATE TABLE `ontology_xref` (
  `object_xref_id` integer  NOT NULL DEFAULT '0'
,  `source_xref_id` integer  DEFAULT NULL
,  `linkage_type` varchar(3) DEFAULT NULL
,  UNIQUE (`object_xref_id`,`source_xref_id`,`linkage_type`)
);
CREATE TABLE `operon` (
  `operon_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `seq_region_id` integer  NOT NULL
,  `seq_region_start` integer  NOT NULL
,  `seq_region_end` integer  NOT NULL
,  `seq_region_strand` integer NOT NULL
,  `display_label` varchar(255) DEFAULT NULL
,  `analysis_id` integer  NOT NULL
,  `stable_id` varchar(128) DEFAULT NULL
,  `version` integer  DEFAULT NULL
,  `created_date` datetime DEFAULT NULL
,  `modified_date` datetime DEFAULT NULL
);
CREATE TABLE `operon_transcript` (
  `operon_transcript_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `seq_region_id` integer  NOT NULL
,  `seq_region_start` integer  NOT NULL
,  `seq_region_end` integer  NOT NULL
,  `seq_region_strand` integer NOT NULL
,  `operon_id` integer  NOT NULL
,  `display_label` varchar(255) DEFAULT NULL
,  `analysis_id` integer  NOT NULL
,  `stable_id` varchar(128) DEFAULT NULL
,  `version` integer  DEFAULT NULL
,  `created_date` datetime DEFAULT NULL
,  `modified_date` datetime DEFAULT NULL
);
CREATE TABLE `operon_transcript_gene` (
  `operon_transcript_id` integer  DEFAULT NULL
,  `gene_id` integer  DEFAULT NULL
);
CREATE TABLE `peptide_archive` (
  `peptide_archive_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `md5_checksum` varchar(32) DEFAULT NULL
,  `peptide_seq` mediumtext NOT NULL
);
CREATE TABLE `prediction_exon` (
  `prediction_exon_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `prediction_transcript_id` integer  NOT NULL
,  `exon_rank` integer  NOT NULL
,  `seq_region_id` integer  NOT NULL
,  `seq_region_start` integer  NOT NULL
,  `seq_region_end` integer  NOT NULL
,  `seq_region_strand` integer NOT NULL
,  `start_phase` integer NOT NULL
,  `score` double DEFAULT NULL
,  `p_value` double DEFAULT NULL
);
CREATE TABLE `prediction_transcript` (
  `prediction_transcript_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `seq_region_id` integer  NOT NULL
,  `seq_region_start` integer  NOT NULL
,  `seq_region_end` integer  NOT NULL
,  `seq_region_strand` integer NOT NULL
,  `analysis_id` integer  NOT NULL
,  `display_label` varchar(255) DEFAULT NULL
);
CREATE TABLE `protein_align_feature` (
  `protein_align_feature_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `seq_region_id` integer  NOT NULL
,  `seq_region_start` integer  NOT NULL
,  `seq_region_end` integer  NOT NULL
,  `seq_region_strand` integer NOT NULL DEFAULT '1'
,  `hit_start` integer NOT NULL
,  `hit_end` integer NOT NULL
,  `hit_name` varchar(40) NOT NULL
,  `analysis_id` integer  NOT NULL
,  `score` double DEFAULT NULL
,  `evalue` double DEFAULT NULL
,  `perc_ident` float DEFAULT NULL
,  `cigar_line` text
,  `external_db_id` integer  DEFAULT NULL
,  `hcoverage` double DEFAULT NULL
,  `align_type` text  DEFAULT 'ensembl'
);
CREATE TABLE `protein_feature` (
  `protein_feature_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `translation_id` integer  NOT NULL
,  `seq_start` integer NOT NULL
,  `seq_end` integer NOT NULL
,  `hit_start` integer NOT NULL
,  `hit_end` integer NOT NULL
,  `hit_name` varchar(40) NOT NULL
,  `analysis_id` integer  NOT NULL
,  `score` double DEFAULT NULL
,  `evalue` double DEFAULT NULL
,  `perc_ident` float DEFAULT NULL
,  `external_data` text
,  `hit_description` text
,  `cigar_line` text
,  `align_type` text  DEFAULT NULL
,  UNIQUE (`translation_id`,`hit_name`,`seq_start`,`seq_end`,`hit_start`,`hit_end`,`analysis_id`)
);
CREATE TABLE `repeat_consensus` (
  `repeat_consensus_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `repeat_name` varchar(255) NOT NULL
,  `repeat_class` varchar(100) NOT NULL
,  `repeat_type` varchar(40) NOT NULL
,  `repeat_consensus` text
);
CREATE TABLE `repeat_feature` (
  `repeat_feature_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `seq_region_id` integer  NOT NULL
,  `seq_region_start` integer  NOT NULL
,  `seq_region_end` integer  NOT NULL
,  `seq_region_strand` integer NOT NULL DEFAULT '1'
,  `repeat_start` integer NOT NULL
,  `repeat_end` integer NOT NULL
,  `repeat_consensus_id` integer  NOT NULL
,  `analysis_id` integer  NOT NULL
,  `score` double DEFAULT NULL
);
CREATE TABLE `rnaproduct` (
  `rnaproduct_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `rnaproduct_type_id` integer  NOT NULL
,  `transcript_id` integer  NOT NULL
,  `seq_start` integer NOT NULL
,  `start_exon_id` integer  DEFAULT NULL
,  `seq_end` integer NOT NULL
,  `end_exon_id` integer  DEFAULT NULL
,  `stable_id` varchar(128) DEFAULT NULL
,  `version` integer  DEFAULT NULL
,  `created_date` datetime DEFAULT NULL
,  `modified_date` datetime DEFAULT NULL
);
CREATE TABLE `rnaproduct_attrib` (
  `rnaproduct_id` integer  NOT NULL DEFAULT '0'
,  `attrib_type_id` integer  NOT NULL DEFAULT '0'
,  `value` text NOT NULL
,  UNIQUE (`rnaproduct_id`,`attrib_type_id`,`value`)
);
CREATE TABLE `rnaproduct_type` (
  `rnaproduct_type_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `code` varchar(20) NOT NULL DEFAULT ''
,  `name` varchar(255) NOT NULL DEFAULT ''
,  `description` text
,  UNIQUE (`code`)
);
CREATE TABLE `seq_region` (
  `seq_region_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `name` varchar(255) NOT NULL
,  `coord_system_id` integer  NOT NULL
,  `length` integer  NOT NULL
,  UNIQUE (`name`,`coord_system_id`)
);
CREATE TABLE `seq_region_attrib` (
  `seq_region_id` integer  NOT NULL DEFAULT '0'
,  `attrib_type_id` integer  NOT NULL DEFAULT '0'
,  `value` text NOT NULL
,  UNIQUE (`seq_region_id`,`attrib_type_id`,`value`)
);
CREATE TABLE `seq_region_mapping` (
  `external_seq_region_id` integer  NOT NULL
,  `internal_seq_region_id` integer  NOT NULL
,  `mapping_set_id` integer  NOT NULL
);
CREATE TABLE `seq_region_synonym` (
  `seq_region_synonym_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `seq_region_id` integer  NOT NULL
,  `synonym` varchar(250) NOT NULL
,  `external_db_id` integer  DEFAULT NULL
,  UNIQUE (`synonym`,`seq_region_id`)
);
CREATE TABLE `simple_feature` (
  `simple_feature_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `seq_region_id` integer  NOT NULL
,  `seq_region_start` integer  NOT NULL
,  `seq_region_end` integer  NOT NULL
,  `seq_region_strand` integer NOT NULL
,  `display_label` varchar(255) NOT NULL
,  `analysis_id` integer  NOT NULL
,  `score` double DEFAULT NULL
);
CREATE TABLE `stable_id_event` (
  `old_stable_id` varchar(128) DEFAULT NULL
,  `old_version` integer DEFAULT NULL
,  `new_stable_id` varchar(128) DEFAULT NULL
,  `new_version` integer DEFAULT NULL
,  `mapping_session_id` integer  NOT NULL DEFAULT '0'
,  `type` text  NOT NULL
,  `score` float NOT NULL DEFAULT '0'
,  UNIQUE (`mapping_session_id`,`old_stable_id`,`new_stable_id`,`type`)
);
CREATE TABLE `supporting_feature` (
  `exon_id` integer  NOT NULL DEFAULT '0'
,  `feature_type` text  DEFAULT NULL
,  `feature_id` integer  NOT NULL DEFAULT '0'
,  UNIQUE (`exon_id`,`feature_type`,`feature_id`)
);
CREATE TABLE `transcript` (
  `transcript_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `gene_id` integer  DEFAULT NULL
,  `analysis_id` integer  NOT NULL
,  `seq_region_id` integer  NOT NULL
,  `seq_region_start` integer  NOT NULL
,  `seq_region_end` integer  NOT NULL
,  `seq_region_strand` integer NOT NULL
,  `display_xref_id` integer  DEFAULT NULL
,  `source` varchar(40) NOT NULL DEFAULT 'ensembl'
,  `biotype` varchar(40) NOT NULL
,  `description` text
,  `is_current` integer NOT NULL DEFAULT '1'
,  `canonical_translation_id` integer  DEFAULT NULL
,  `stable_id` varchar(128) DEFAULT NULL
,  `version` integer  DEFAULT NULL
,  `created_date` datetime DEFAULT NULL
,  `modified_date` datetime DEFAULT NULL
,  UNIQUE (`canonical_translation_id`)
);
CREATE TABLE `transcript_attrib` (
  `transcript_id` integer  NOT NULL DEFAULT '0'
,  `attrib_type_id` integer  NOT NULL DEFAULT '0'
,  `value` text NOT NULL
,  UNIQUE (`transcript_id`,`attrib_type_id`,`value`)
);
CREATE TABLE `transcript_intron_supporting_evidence` (
  `transcript_id` integer  NOT NULL
,  `intron_supporting_evidence_id` integer  NOT NULL
,  `previous_exon_id` integer  NOT NULL
,  `next_exon_id` integer  NOT NULL
,  PRIMARY KEY (`intron_supporting_evidence_id`,`transcript_id`)
);
CREATE TABLE `transcript_supporting_feature` (
  `transcript_id` integer  NOT NULL DEFAULT '0'
,  `feature_type` text  DEFAULT NULL
,  `feature_id` integer  NOT NULL DEFAULT '0'
,  UNIQUE (`transcript_id`,`feature_type`,`feature_id`)
);
CREATE TABLE `translation` (
  `translation_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `transcript_id` integer  NOT NULL
,  `seq_start` integer NOT NULL
,  `start_exon_id` integer  NOT NULL
,  `seq_end` integer NOT NULL
,  `end_exon_id` integer  NOT NULL
,  `stable_id` varchar(128) DEFAULT NULL
,  `version` integer  DEFAULT NULL
,  `created_date` datetime DEFAULT NULL
,  `modified_date` datetime DEFAULT NULL
);
CREATE TABLE `translation_attrib` (
  `translation_id` integer  NOT NULL DEFAULT '0'
,  `attrib_type_id` integer  NOT NULL DEFAULT '0'
,  `value` text NOT NULL
,  UNIQUE (`translation_id`,`attrib_type_id`,`value`)
);
CREATE TABLE `unmapped_object` (
  `unmapped_object_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `type` text  NOT NULL
,  `analysis_id` integer  NOT NULL
,  `external_db_id` integer  DEFAULT NULL
,  `identifier` varchar(255) NOT NULL
,  `unmapped_reason_id` integer  NOT NULL
,  `query_score` double DEFAULT NULL
,  `target_score` double DEFAULT NULL
,  `ensembl_id` integer  DEFAULT '0'
,  `ensembl_object_type` text  DEFAULT 'RawContig'
,  `parent` varchar(255) DEFAULT NULL
,  UNIQUE (`ensembl_id`,`ensembl_object_type`,`identifier`,`unmapped_reason_id`,`parent`,`external_db_id`)
);
CREATE TABLE `unmapped_reason` (
  `unmapped_reason_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `summary_description` varchar(255) DEFAULT NULL
,  `full_description` varchar(255) DEFAULT NULL
);
CREATE TABLE `xref` (
  `xref_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `external_db_id` integer  NOT NULL
,  `dbprimary_acc` varchar(512) NOT NULL
,  `display_label` varchar(512) NOT NULL
,  `version` varchar(10) DEFAULT NULL
,  `description` text
,  `info_type` text  NOT NULL DEFAULT 'NONE'
,  `info_text` varchar(255) NOT NULL DEFAULT ''
,  UNIQUE (`dbprimary_acc`,`external_db_id`,`info_type`,`info_text`,`version`)
);
CREATE INDEX "idx_translation_attrib_type_val_idx" ON "translation_attrib" (`attrib_type_id`,`value`);
CREATE INDEX "idx_translation_attrib_val_only_idx" ON "translation_attrib" (`value`);
CREATE INDEX "idx_translation_attrib_translation_idx" ON "translation_attrib" (`translation_id`);
CREATE INDEX "idx_prediction_exon_transcript_idx" ON "prediction_exon" (`prediction_transcript_id`);
CREATE INDEX "idx_prediction_exon_seq_region_idx" ON "prediction_exon" (`seq_region_id`,`seq_region_start`);
CREATE INDEX "idx_operon_transcript_operon_idx" ON "operon_transcript" (`operon_id`);
CREATE INDEX "idx_operon_transcript_seq_region_idx" ON "operon_transcript" (`seq_region_id`,`seq_region_start`);
CREATE INDEX "idx_operon_transcript_stable_id_idx" ON "operon_transcript" (`stable_id`,`version`);
CREATE INDEX "idx_marker_synonym_marker_synonym_idx" ON "marker_synonym" (`marker_synonym_id`,`name`);
CREATE INDEX "idx_marker_synonym_marker_idx" ON "marker_synonym" (`marker_id`);
CREATE INDEX "idx_dna_align_feature_attrib_dna_align_feature_idx" ON "dna_align_feature_attrib" (`dna_align_feature_id`);
CREATE INDEX "idx_dna_align_feature_attrib_type_val_idx" ON "dna_align_feature_attrib" (`attrib_type_id`,`value`);
CREATE INDEX "idx_dna_align_feature_attrib_val_only_idx" ON "dna_align_feature_attrib" (`value`);
CREATE INDEX "idx_data_file_df_name_idx" ON "data_file" (`name`);
CREATE INDEX "idx_data_file_df_analysis_idx" ON "data_file" (`analysis_id`);
CREATE INDEX "idx_misc_feature_seq_region_idx" ON "misc_feature" (`seq_region_id`,`seq_region_start`);
CREATE INDEX "idx_protein_align_feature_seq_region_idx" ON "protein_align_feature" (`seq_region_id`,`analysis_id`,`seq_region_start`,`score`);
CREATE INDEX "idx_protein_align_feature_seq_region_idx_2" ON "protein_align_feature" (`seq_region_id`,`seq_region_start`);
CREATE INDEX "idx_protein_align_feature_hit_idx" ON "protein_align_feature" (`hit_name`);
CREATE INDEX "idx_protein_align_feature_analysis_idx" ON "protein_align_feature" (`analysis_id`);
CREATE INDEX "idx_protein_align_feature_external_db_idx" ON "protein_align_feature" (`external_db_id`);
CREATE INDEX "idx_associated_xref_associated_source_idx" ON "associated_xref" (`source_xref_id`);
CREATE INDEX "idx_associated_xref_associated_object_idx" ON "associated_xref" (`object_xref_id`);
CREATE INDEX "idx_associated_xref_associated_idx" ON "associated_xref" (`xref_id`);
CREATE INDEX "idx_associated_xref_associated_group_idx" ON "associated_xref" (`associated_group_id`);
CREATE INDEX "idx_xref_display_index" ON "xref" (`display_label`);
CREATE INDEX "idx_xref_info_type_idx" ON "xref" (`info_type`);
CREATE INDEX "idx_transcript_intron_supporting_evidence_transcript_idx" ON "transcript_intron_supporting_evidence" (`transcript_id`);
CREATE INDEX "idx_operon_transcript_gene_operon_transcript_gene_idx" ON "operon_transcript_gene" (`operon_transcript_id`,`gene_id`);
CREATE INDEX "idx_seq_region_attrib_type_val_idx" ON "seq_region_attrib" (`attrib_type_id`,`value`);
CREATE INDEX "idx_seq_region_attrib_val_only_idx" ON "seq_region_attrib" (`value`);
CREATE INDEX "idx_seq_region_attrib_seq_region_idx" ON "seq_region_attrib" (`seq_region_id`);
CREATE INDEX "idx_operon_seq_region_idx" ON "operon" (`seq_region_id`,`seq_region_start`);
CREATE INDEX "idx_operon_name_idx" ON "operon" (`display_label`);
CREATE INDEX "idx_operon_stable_id_idx" ON "operon" (`stable_id`,`version`);
CREATE INDEX "idx_stable_id_event_new_idx" ON "stable_id_event" (`new_stable_id`);
CREATE INDEX "idx_stable_id_event_old_idx" ON "stable_id_event" (`old_stable_id`);
CREATE INDEX "idx_repeat_feature_seq_region_idx" ON "repeat_feature" (`seq_region_id`,`seq_region_start`);
CREATE INDEX "idx_repeat_feature_repeat_idx" ON "repeat_feature" (`repeat_consensus_id`);
CREATE INDEX "idx_repeat_feature_analysis_idx" ON "repeat_feature" (`analysis_id`);
CREATE INDEX "idx_dependent_xref_dependent" ON "dependent_xref" (`dependent_xref_id`);
CREATE INDEX "idx_dependent_xref_master_idx" ON "dependent_xref" (`master_xref_id`);
CREATE INDEX "idx_density_feature_seq_region_idx" ON "density_feature" (`density_type_id`,`seq_region_id`,`seq_region_start`);
CREATE INDEX "idx_density_feature_seq_region_id_idx" ON "density_feature" (`seq_region_id`);
CREATE INDEX "idx_peptide_archive_checksum" ON "peptide_archive" (`md5_checksum`);
CREATE INDEX "idx_intron_supporting_evidence_seq_region_idx" ON "intron_supporting_evidence" (`seq_region_id`,`seq_region_start`);
CREATE INDEX "idx_exon_seq_region_idx" ON "exon" (`seq_region_id`,`seq_region_start`);
CREATE INDEX "idx_exon_stable_id_idx" ON "exon" (`stable_id`,`version`);
CREATE INDEX "idx_ontology_xref_source_idx" ON "ontology_xref" (`source_xref_id`);
CREATE INDEX "idx_ontology_xref_object_idx" ON "ontology_xref" (`object_xref_id`);
CREATE INDEX "idx_assembly_cmp_seq_region_idx" ON "assembly" (`cmp_seq_region_id`);
CREATE INDEX "idx_assembly_asm_seq_region_idx" ON "assembly" (`asm_seq_region_id`,`asm_start`);
CREATE INDEX "idx_exon_transcript_transcript" ON "exon_transcript" (`transcript_id`);
CREATE INDEX "idx_exon_transcript_exon" ON "exon_transcript" (`exon_id`);
CREATE INDEX "idx_alt_allele_attrib_aa_idx" ON "alt_allele_attrib" (`alt_allele_id`,`attrib`);
CREATE INDEX "idx_unmapped_object_id_idx" ON "unmapped_object" (`identifier`);
CREATE INDEX "idx_unmapped_object_anal_exdb_idx" ON "unmapped_object" (`analysis_id`,`external_db_id`);
CREATE INDEX "idx_unmapped_object_ext_db_identifier_idx" ON "unmapped_object" (`external_db_id`,`identifier`);
CREATE INDEX "idx_seq_region_synonym_seq_region_idx" ON "seq_region_synonym" (`seq_region_id`);
CREATE INDEX "idx_rnaproduct_attrib_type_val_idx" ON "rnaproduct_attrib" (`attrib_type_id`,`value`);
CREATE INDEX "idx_rnaproduct_attrib_val_only_idx" ON "rnaproduct_attrib" (`value`);
CREATE INDEX "idx_rnaproduct_attrib_rnaproduct_idx" ON "rnaproduct_attrib" (`rnaproduct_id`);
CREATE INDEX "idx_object_xref_ensembl_idx" ON "object_xref" (`ensembl_object_type`,`ensembl_id`);
CREATE INDEX "idx_object_xref_analysis_idx" ON "object_xref" (`analysis_id`);
CREATE INDEX "idx_ditag_feature_ditag_idx" ON "ditag_feature" (`ditag_id`);
CREATE INDEX "idx_ditag_feature_ditag_pair_idx" ON "ditag_feature" (`ditag_pair_id`);
CREATE INDEX "idx_ditag_feature_seq_region_idx" ON "ditag_feature" (`seq_region_id`,`seq_region_start`,`seq_region_end`);
CREATE INDEX "idx_repeat_consensus_name" ON "repeat_consensus" (`repeat_name`);
CREATE INDEX "idx_repeat_consensus_class" ON "repeat_consensus" (`repeat_class`);
CREATE INDEX "idx_repeat_consensus_consensus" ON "repeat_consensus" (`repeat_consensus`);
CREATE INDEX "idx_repeat_consensus_type" ON "repeat_consensus" (`repeat_type`);
CREATE INDEX "idx_marker_map_location_map_idx" ON "marker_map_location" (`map_id`,`chromosome_name`,`position`);
CREATE INDEX "idx_karyotype_region_band_idx" ON "karyotype" (`seq_region_id`,`band`);
CREATE INDEX "idx_meta_species_value_idx" ON "meta" (`species_id`,`meta_value`);
CREATE INDEX "idx_protein_feature_translation_idx" ON "protein_feature" (`translation_id`);
CREATE INDEX "idx_protein_feature_hitname_idx" ON "protein_feature" (`hit_name`);
CREATE INDEX "idx_protein_feature_analysis_idx" ON "protein_feature" (`analysis_id`);
CREATE INDEX "idx_gene_archive_gene_idx" ON "gene_archive" (`gene_stable_id`,`gene_version`);
CREATE INDEX "idx_gene_archive_transcript_idx" ON "gene_archive" (`transcript_stable_id`,`transcript_version`);
CREATE INDEX "idx_gene_archive_translation_idx" ON "gene_archive" (`translation_stable_id`,`translation_version`);
CREATE INDEX "idx_gene_archive_peptide_archive_id_idx" ON "gene_archive" (`peptide_archive_id`);
CREATE INDEX "idx_dna_align_feature_seq_region_idx" ON "dna_align_feature" (`seq_region_id`,`analysis_id`,`seq_region_start`,`score`);
CREATE INDEX "idx_dna_align_feature_seq_region_idx_2" ON "dna_align_feature" (`seq_region_id`,`seq_region_start`);
CREATE INDEX "idx_dna_align_feature_hit_idx" ON "dna_align_feature" (`hit_name`);
CREATE INDEX "idx_dna_align_feature_analysis_idx" ON "dna_align_feature" (`analysis_id`);
CREATE INDEX "idx_dna_align_feature_external_db_idx" ON "dna_align_feature" (`external_db_id`);
CREATE INDEX "idx_rnaproduct_transcript_idx" ON "rnaproduct" (`transcript_id`);
CREATE INDEX "idx_rnaproduct_stable_id_idx" ON "rnaproduct" (`stable_id`,`version`);
CREATE INDEX "idx_transcript_supporting_feature_feature_idx" ON "transcript_supporting_feature" (`feature_type`,`feature_id`);
CREATE INDEX "idx_assembly_exception_sr_idx" ON "assembly_exception" (`seq_region_id`,`seq_region_start`);
CREATE INDEX "idx_assembly_exception_ex_idx" ON "assembly_exception" (`exc_seq_region_id`,`exc_seq_region_start`);
CREATE INDEX "idx_supporting_feature_feature_idx" ON "supporting_feature" (`feature_type`,`feature_id`);
CREATE INDEX "idx_gene_attrib_type_val_idx" ON "gene_attrib" (`attrib_type_id`,`value`);
CREATE INDEX "idx_gene_attrib_val_only_idx" ON "gene_attrib" (`value`);
CREATE INDEX "idx_gene_attrib_gene_idx" ON "gene_attrib" (`gene_id`);
CREATE INDEX "idx_coord_system_species_idx" ON "coord_system" (`species_id`);
CREATE INDEX "idx_seq_region_mapping_mapping_set_idx" ON "seq_region_mapping" (`mapping_set_id`);
CREATE INDEX "idx_misc_feature_misc_set_reverse_idx" ON "misc_feature_misc_set" (`misc_set_id`,`misc_feature_id`);
CREATE INDEX "idx_interpro_id_idx" ON "interpro" (`id`);
CREATE INDEX "idx_alt_allele_gene_id" ON "alt_allele" (`gene_id`,`alt_allele_group_id`);
CREATE INDEX "idx_transcript_seq_region_idx" ON "transcript" (`seq_region_id`,`seq_region_start`);
CREATE INDEX "idx_transcript_gene_index" ON "transcript" (`gene_id`);
CREATE INDEX "idx_transcript_xref_id_index" ON "transcript" (`display_xref_id`);
CREATE INDEX "idx_transcript_analysis_idx" ON "transcript" (`analysis_id`);
CREATE INDEX "idx_transcript_stable_id_idx" ON "transcript" (`stable_id`,`version`);
CREATE INDEX "idx_transcript_attrib_type_val_idx" ON "transcript_attrib" (`attrib_type_id`,`value`);
CREATE INDEX "idx_transcript_attrib_val_only_idx" ON "transcript_attrib" (`value`);
CREATE INDEX "idx_transcript_attrib_transcript_idx" ON "transcript_attrib" (`transcript_id`);
CREATE INDEX "idx_prediction_transcript_seq_region_idx" ON "prediction_transcript" (`seq_region_id`,`seq_region_start`);
CREATE INDEX "idx_prediction_transcript_analysis_idx" ON "prediction_transcript" (`analysis_id`);
CREATE INDEX "idx_marker_marker_idx" ON "marker" (`marker_id`,`priority`);
CREATE INDEX "idx_marker_display_idx" ON "marker" (`display_marker_synonym_id`);
CREATE INDEX "idx_translation_transcript_idx" ON "translation" (`transcript_id`);
CREATE INDEX "idx_translation_stable_id_idx" ON "translation" (`stable_id`,`version`);
CREATE INDEX "idx_simple_feature_seq_region_idx" ON "simple_feature" (`seq_region_id`,`seq_region_start`);
CREATE INDEX "idx_simple_feature_analysis_idx" ON "simple_feature" (`analysis_id`);
CREATE INDEX "idx_simple_feature_hit_idx" ON "simple_feature" (`display_label`);
CREATE INDEX "idx_gene_seq_region_idx" ON "gene" (`seq_region_id`,`seq_region_start`);
CREATE INDEX "idx_gene_xref_id_index" ON "gene" (`display_xref_id`);
CREATE INDEX "idx_gene_analysis_idx" ON "gene" (`analysis_id`);
CREATE INDEX "idx_gene_stable_id_idx" ON "gene" (`stable_id`,`version`);
CREATE INDEX "idx_gene_canonical_transcript_id_idx" ON "gene" (`canonical_transcript_id`);
CREATE INDEX "idx_external_synonym_name_index" ON "external_synonym" (`synonym`);
CREATE INDEX "idx_seq_region_cs_idx" ON "seq_region" (`coord_system_id`);
CREATE INDEX "idx_misc_attrib_type_val_idx" ON "misc_attrib" (`attrib_type_id`,`value`);
CREATE INDEX "idx_misc_attrib_val_only_idx" ON "misc_attrib" (`value`);
CREATE INDEX "idx_misc_attrib_misc_feature_idx" ON "misc_attrib" (`misc_feature_id`);
CREATE INDEX "idx_marker_feature_seq_region_idx" ON "marker_feature" (`seq_region_id`,`seq_region_start`);
CREATE INDEX "idx_marker_feature_analysis_idx" ON "marker_feature" (`analysis_id`);
