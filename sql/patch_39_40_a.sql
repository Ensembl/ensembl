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

# patch_39_40_a
#
# title: rationalise key columns
#;
# description:
# rationalise all the primary and foreign key columns to be INT(10) UNSIGNED
# Note that this is not a definitive list of all primary/foreign keys, just
# the ones that were not defined as INT(10) UNSIGNED in schema 39.

ALTER TABLE oligo_feature CHANGE COLUMN oligo_feature_id oligo_feature_id INT(10) UNSIGNED NOT NULL AUTO_INCREMENT;
ALTER TABLE oligo_feature CHANGE COLUMN seq_region_id seq_region_id INT(10) UNSIGNED NOT NULL;
ALTER TABLE oligo_feature CHANGE COLUMN oligo_probe_id oligo_probe_id INT(10) UNSIGNED NOT NULL;

ALTER TABLE oligo_probe CHANGE COLUMN oligo_probe_id oligo_probe_id INT(10) UNSIGNED NOT NULL AUTO_INCREMENT;
ALTER TABLE oligo_probe CHANGE COLUMN oligo_array_id oligo_array_id INT(10) UNSIGNED NOT NULL;

ALTER TABLE oligo_array CHANGE COLUMN oligo_array_id oligo_array_id INT(10) UNSIGNED NOT NULL auto_increment;
ALTER TABLE oligo_array CHANGE COLUMN parent_array_id parent_array_id INT(10) UNSIGNED;

ALTER TABLE alt_allele CHANGE COLUMN alt_allele_id alt_allele_id INT(10) UNSIGNED NOT NULL AUTO_INCREMENT;
ALTER TABLE alt_allele CHANGE COLUMN gene_id gene_id INT(10) UNSIGNED NOT NULL;

ALTER TABLE dna CHANGE COLUMN seq_region_id seq_region_id INT(10) UNSIGNED NOT NULL;

ALTER TABLE dnac CHANGE COLUMN seq_region_id seq_region_id INT(10) UNSIGNED NOT NULL;

ALTER TABLE exon CHANGE COLUMN exon_id exon_id INT(10) UNSIGNED NOT NULL AUTO_INCREMENT;

ALTER TABLE exon_stable_id CHANGE COLUMN exon_id exon_id INT(10) UNSIGNED NOT NULL;

ALTER TABLE exon_transcript CHANGE COLUMN exon_id exon_id INT(10) UNSIGNED NOT NULL;
ALTER TABLE exon_transcript CHANGE COLUMN transcript_id transcript_id INT(10) UNSIGNED NOT NULL;

ALTER TABLE simple_feature CHANGE COLUMN simple_feature_id simple_feature_id INT(10) UNSIGNED NOT NULL AUTO_INCREMENT;

ALTER TABLE protein_align_feature CHANGE COLUMN protein_align_feature_id protein_align_feature_id INT(10) UNSIGNED NOT NULL AUTO_INCREMENT;

ALTER TABLE dna_align_feature CHANGE COLUMN dna_align_feature_id dna_align_feature_id INT(10) UNSIGNED NOT NULL AUTO_INCREMENT;

ALTER TABLE repeat_consensus CHANGE COLUMN repeat_consensus_id repeat_consensus_id INT(10) UNSIGNED NOT NULL AUTO_INCREMENT;

ALTER TABLE repeat_feature CHANGE COLUMN repeat_feature_id repeat_feature_id INT(10) UNSIGNED NOT NULL AUTO_INCREMENT;

ALTER TABLE gene CHANGE COLUMN gene_id gene_id INT(10) UNSIGNED NOT NULL AUTO_INCREMENT;
ALTER TABLE gene CHANGE COLUMN display_xref_id display_xref_id INT(10) UNSIGNED;

ALTER TABLE gene_stable_id CHANGE COLUMN gene_id gene_id INT UNSIGNED NOT NULL;

ALTER TABLE supporting_feature CHANGE COLUMN exon_id exon_id INT(10) UNSIGNED DEFAULT '0' NOT NULL;
ALTER TABLE supporting_feature CHANGE COLUMN feature_id feature_id INT(10) UNSIGNED DEFAULT '0' NOT NULL;

ALTER TABLE transcript_supporting_feature CHANGE COLUMN transcript_id transcript_id INT(10) UNSIGNED DEFAULT '0' NOT NULL;
ALTER TABLE transcript_supporting_feature CHANGE COLUMN feature_id feature_id INT(10) UNSIGNED DEFAULT '0' NOT NULL;

ALTER TABLE transcript CHANGE COLUMN transcript_id transcript_id INT(10) UNSIGNED NOT NULL AUTO_INCREMENT;
ALTER TABLE transcript CHANGE COLUMN gene_id gene_id INT(10) UNSIGNED NOT NULL;
ALTER TABLE transcript CHANGE COLUMN display_xref_id display_xref_id INT(10) UNSIGNED;

ALTER TABLE transcript_stable_id CHANGE COLUMN transcript_id transcript_id INT(10) UNSIGNED NOT NULL;

ALTER TABLE translation CHANGE COLUMN translation_id translation_id INT(10) UNSIGNED NOT NULL AUTO_INCREMENT ;
ALTER TABLE translation CHANGE COLUMN transcript_id transcript_id INT(10) UNSIGNED NOT NULL;
ALTER TABLE translation CHANGE COLUMN start_exon_id start_exon_id INT(10) UNSIGNED NOT NULL;
ALTER TABLE translation CHANGE COLUMN end_exon_id end_exon_id INT(10) UNSIGNED NOT NULL;

ALTER TABLE translation_stable_id CHANGE COLUMN translation_id translation_id INT(10) UNSIGNED NOT NULL;

ALTER TABLE assembly CHANGE COLUMN asm_seq_region_id asm_seq_region_id INT(10) UNSIGNED NOT NULL;

ALTER TABLE protein_feature CHANGE COLUMN translation_id translation_id INT(10) UNSIGNED NOT NULL;

ALTER TABLE karyotype CHANGE COLUMN karyotype_id karyotype_id INT(10) UNSIGNED NOT NULL AUTO_INCREMENT;
ALTER TABLE karyotype CHANGE COLUMN seq_region_id seq_region_id INT(10) UNSIGNED NOT NULL;

ALTER TABLE object_xref CHANGE COLUMN object_xref_id object_xref_id INT(10) UNSIGNED NOT NULL AUTO_INCREMENT;
ALTER TABLE object_xref CHANGE COLUMN ensembl_id ensembl_id INT(10) UNSIGNED NOT NULL;

ALTER TABLE identity_xref CHANGE COLUMN object_xref_id object_xref_id INT(10) UNSIGNED NOT NULL;

ALTER TABLE xref CHANGE COLUMN xref_id xref_id INT(10) UNSIGNED NOT NULL AUTO_INCREMENT;
ALTER TABLE xref CHANGE COLUMN external_db_id external_db_id INT(10) UNSIGNED NOT NULL;

ALTER TABLE external_synonym CHANGE COLUMN xref_id xref_id INT(10) UNSIGNED NOT NULL;

ALTER TABLE external_db CHANGE COLUMN external_db_id external_db_id	 INT(10) UNSIGNED NOT NULL;

ALTER TABLE prediction_exon CHANGE COLUMN prediction_exon_id prediction_exon_id INT(10) UNSIGNED NOT NULL AUTO_INCREMENT;
ALTER TABLE prediction_exon CHANGE COLUMN prediction_transcript_id prediction_transcript_id INT(10) UNSIGNED NOT NULL;
ALTER TABLE prediction_exon CHANGE COLUMN seq_region_id seq_region_id INT(10) UNSIGNED NOT NULL;

ALTER TABLE prediction_transcript CHANGE COLUMN prediction_transcript_id prediction_transcript_id INT(10) UNSIGNED NOT NULL AUTO_INCREMENT;
ALTER TABLE prediction_transcript CHANGE COLUMN seq_region_id seq_region_id INT(10) UNSIGNED NOT NULL;

ALTER TABLE marker_synonym CHANGE COLUMN marker_synonym_id marker_synonym_id INT(10) UNSIGNED NOT NULL AUTO_INCREMENT;
ALTER TABLE marker_synonym CHANGE COLUMN marker_id marker_id INT(10) UNSIGNED NOT NULL;

ALTER TABLE marker CHANGE COLUMN marker_id marker_id INT(10) UNSIGNED NOT NULL AUTO_INCREMENT;
ALTER TABLE marker CHANGE COLUMN display_marker_synonym_id display_marker_synonym_id INT(10) UNSIGNED;

ALTER TABLE marker_feature CHANGE COLUMN marker_feature_id marker_feature_id INT(10) UNSIGNED NOT NULL AUTO_INCREMENT;
ALTER TABLE marker_feature CHANGE COLUMN marker_id marker_id INT(10) UNSIGNED NOT NULL;

ALTER TABLE marker_map_location CHANGE COLUMN marker_id marker_id INT UNSIGNED NOT NULL;
ALTER TABLE marker_map_location CHANGE COLUMN map_id map_id INT UNSIGNED NOT NULL;

ALTER TABLE marker_map_location CHANGE COLUMN marker_id marker_id INT(10) UNSIGNED NOT NULL;
ALTER TABLE marker_map_location CHANGE COLUMN map_id map_id INT(10) UNSIGNED NOT NULL;
ALTER TABLE marker_map_location CHANGE COLUMN marker_synonym_id marker_synonym_id INT(10) UNSIGNED NOT NULL;

ALTER TABLE map CHANGE COLUMN map_id map_id INT(10) UNSIGNED NOT NULL AUTO_INCREMENT;

ALTER TABLE qtl CHANGE COLUMN qtl_id qtl_id INT(10) UNSIGNED AUTO_INCREMENT NOT NULL;
ALTER TABLE qtl CHANGE COLUMN flank_marker_id_1 flank_marker_id_1 INT(10) UNSIGNED;
ALTER TABLE qtl CHANGE COLUMN flank_marker_id_2 flank_marker_id_2 INT(10) UNSIGNED;
ALTER TABLE qtl CHANGE COLUMN peak_marker_id peak_marker_id INT(10) UNSIGNED;

ALTER TABLE qtl_synonym CHANGE COLUMN qtl_synonym_id qtl_synonym_id 	 INT(10) UNSIGNED AUTO_INCREMENT NOT NULL;
ALTER TABLE qtl_synonym CHANGE COLUMN qtl_id qtl_id INT(10) UNSIGNED NOT NULL;

ALTER TABLE qtl_feature CHANGE COLUMN seq_region_id seq_region_id 	INT(10)	UNSIGNED NOT NULL;
ALTER TABLE qtl_feature CHANGE COLUMN qtl_id qtl_id INT(10)	UNSIGNED NOT NULL;

ALTER TABLE mapping_session CHANGE COLUMN mapping_session_id mapping_session_id 	 INT(10) UNSIGNED NOT NULL AUTO_INCREMENT;

ALTER TABLE stable_id_event CHANGE COLUMN mapping_session_id mapping_session_id INT(10) UNSIGNED NOT NULL DEFAULT '0';

ALTER TABLE gene_archive CHANGE COLUMN peptide_archive_id peptide_archive_id INT(10) UNSIGNED NOT NULL;
ALTER TABLE gene_archive CHANGE COLUMN mapping_session_id mapping_session_id INT(10) UNSIGNED NOT NULL;

ALTER TABLE peptide_archive CHANGE COLUMN peptide_archive_id peptide_archive_id INT(10) UNSIGNED NOT NULL AUTO_INCREMENT;

ALTER TABLE seq_region CHANGE COLUMN coord_system_id coord_system_id INT(10) UNSIGNED NOT NULL;

ALTER TABLE assembly_exception CHANGE COLUMN seq_region_id seq_region_id INT(10) UNSIGNED NOT NULL;
ALTER TABLE assembly_exception CHANGE COLUMN exc_seq_region_id exc_seq_region_id INT(10) UNSIGNED NOT NULL;

ALTER TABLE coord_system CHANGE COLUMN coord_system_id coord_system_id INT(10) UNSIGNED NOT NULL AUTO_INCREMENT;

ALTER TABLE meta_coord CHANGE COLUMN coord_system_id coord_system_id INT(10) UNSIGNED NOT NULL;

ALTER TABLE density_feature CHANGE COLUMN density_feature_id density_feature_id INT(10) UNSIGNED NOT NULL AUTO_INCREMENT;
ALTER TABLE density_feature CHANGE COLUMN density_type_id density_type_id INT(10) UNSIGNED NOT NULL;
ALTER TABLE density_feature CHANGE COLUMN seq_region_id seq_region_id INT(10) UNSIGNED NOT NULL;
ALTER TABLE density_feature CHANGE COLUMN seq_region_start seq_region_start INT(10) UNSIGNED NOT NULL;
ALTER TABLE density_feature CHANGE COLUMN seq_region_end seq_region_end INT(10) UNSIGNED NOT NULL;

ALTER TABLE density_type CHANGE COLUMN density_type_id density_type_id INT(10) UNSIGNED NOT NULL AUTO_INCREMENT;

ALTER TABLE regulatory_feature CHANGE COLUMN regulatory_feature_id regulatory_feature_id INT(10) UNSIGNED NOT NULL AUTO_INCREMENT;
ALTER TABLE regulatory_feature CHANGE COLUMN seq_region_id seq_region_id INT(10) UNSIGNED NOT NULL;
ALTER TABLE regulatory_feature CHANGE COLUMN regulatory_factor_id regulatory_factor_id INT(10) UNSIGNED;

ALTER TABLE regulatory_factor CHANGE COLUMN regulatory_factor_id regulatory_factor_id INT(10) UNSIGNED NOT NULL AUTO_INCREMENT ;

ALTER TABLE regulatory_feature_object CHANGE COLUMN regulatory_feature_id regulatory_feature_id INT(10) UNSIGNED NOT NULL;
ALTER TABLE regulatory_feature_object CHANGE COLUMN ensembl_object_id ensembl_object_id INT(10) UNSIGNED NOT NULL;

ALTER TABLE regulatory_factor_coding CHANGE COLUMN regulatory_factor_id regulatory_factor_id INT(10) UNSIGNED NOT NULL;
ALTER TABLE regulatory_factor_coding CHANGE COLUMN transcript_id transcript_id INT(10) UNSIGNED;
ALTER TABLE regulatory_factor_coding CHANGE COLUMN gene_id gene_id 	INT(10) UNSIGNED;

ALTER TABLE regulatory_search_region CHANGE COLUMN regulatory_search_region_id regulatory_search_region_id INT(10) UNSIGNED NOT NULL AUTO_INCREMENT;
ALTER TABLE regulatory_search_region CHANGE COLUMN seq_region_id seq_region_id INT(10) UNSIGNED NOT NULL;
ALTER TABLE regulatory_search_region CHANGE COLUMN ensembl_object_id ensembl_object_id INT(10) UNSIGNED;

ALTER TABLE unmapped_object CHANGE COLUMN unmapped_object_id unmapped_object_id INT(10) UNSIGNED NOT NULL AUTO_INCREMENT;
ALTER TABLE unmapped_object CHANGE COLUMN external_db_id external_db_id INT(10) UNSIGNED NOT NULL;

# Similarly make all seq_region_start & seq_region_end columns INT(10) UNSIGNED
# Although these are not used as keys, they are used in many joins and having
# the same column type should make joins faster as MySQL will not have to do
# any type casting.

ALTER TABLE oligo_feature CHANGE COLUMN seq_region_start seq_region_start INT(10) UNSIGNED NOT NULL;
ALTER TABLE oligo_feature CHANGE COLUMN seq_region_end seq_region_end INT(10) UNSIGNED NOT NULL;

ALTER TABLE prediction_exon CHANGE COLUMN seq_region_start seq_region_start INT(10) UNSIGNED NOT NULL;
ALTER TABLE prediction_exon CHANGE COLUMN seq_region_end seq_region_end INT(10) UNSIGNED NOT NULL;

ALTER TABLE qtl_feature CHANGE COLUMN seq_region_start seq_region_start INT(10)	UNSIGNED NOT NULL;
ALTER TABLE qtl_feature CHANGE COLUMN seq_region_end seq_region_end INT(10)	UNSIGNED NOT NULL;

ALTER TABLE assembly_exception CHANGE COLUMN seq_region_start seq_region_start INT(10) UNSIGNED NOT NULL;
ALTER TABLE assembly_exception CHANGE COLUMN seq_region_end seq_region_end INT(10) UNSIGNED NOT NULL;
ALTER TABLE assembly_exception CHANGE COLUMN exc_seq_region_start exc_seq_region_start INT(10) UNSIGNED NOT NULL;
ALTER TABLE assembly_exception CHANGE COLUMN exc_seq_region_end exc_seq_region_end INT(10) UNSIGNED NOT NULL;

ALTER TABLE regulatory_feature CHANGE COLUMN seq_region_start seq_region_start INT(10) UNSIGNED NOT NULL;
ALTER TABLE regulatory_feature CHANGE COLUMN seq_region_end seq_region_end INT(10) UNSIGNED NOT NULL;

ALTER TABLE regulatory_search_region CHANGE COLUMN seq_region_start seq_region_start INT(10) UNSIGNED NOT NULL;
ALTER TABLE regulatory_search_region CHANGE COLUMN seq_region_end seq_region_end INT(10) UNSIGNED NOT NULL;

ALTER TABLE prediction_transcript CHANGE COLUMN seq_region_start seq_region_start INT(10) UNSIGNED NOT NULL;
ALTER TABLE prediction_transcript CHANGE COLUMN seq_region_end seq_region_end INT(10) UNSIGNED NOT NULL;

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_39_40_a.sql|rationalise_key_columns');

