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

# patch_65_66_f.sql

# Title: Removal of default options from a number of tables
#
# Description:
# A change in healthchecks has meant the checking of default values which was
# originally changed back in 2006. Programatically generated from the results
# of CompareSchema healthcheck for canis e!66 vs. master_schema_66

ALTER TABLE alt_allele 
   ALTER COLUMN gene_id DROP DEFAULT;

ALTER TABLE analysis_description 
   ALTER COLUMN analysis_id DROP DEFAULT;

ALTER TABLE assembly 
   ALTER COLUMN asm_end DROP DEFAULT,
   ALTER COLUMN asm_seq_region_id DROP DEFAULT,
   ALTER COLUMN asm_start DROP DEFAULT,
   ALTER COLUMN cmp_end DROP DEFAULT,
   ALTER COLUMN cmp_seq_region_id DROP DEFAULT,
   ALTER COLUMN cmp_start DROP DEFAULT,
   ALTER COLUMN ori DROP DEFAULT;

ALTER TABLE assembly_exception 
   ALTER COLUMN exc_seq_region_end DROP DEFAULT,
   ALTER COLUMN exc_seq_region_id DROP DEFAULT,
   ALTER COLUMN exc_seq_region_start DROP DEFAULT,
   ALTER COLUMN ori DROP DEFAULT,
   ALTER COLUMN seq_region_end DROP DEFAULT,
   ALTER COLUMN seq_region_id DROP DEFAULT,
   ALTER COLUMN seq_region_start DROP DEFAULT;

ALTER TABLE coord_system 
   ALTER COLUMN name DROP DEFAULT,
   ALTER COLUMN rank DROP DEFAULT;

ALTER TABLE density_feature 
   ALTER COLUMN density_type_id DROP DEFAULT,
   ALTER COLUMN density_value DROP DEFAULT,
   ALTER COLUMN seq_region_end DROP DEFAULT,
   ALTER COLUMN seq_region_id DROP DEFAULT,
   ALTER COLUMN seq_region_start DROP DEFAULT;

ALTER TABLE density_type 
   ALTER COLUMN analysis_id DROP DEFAULT,
   ALTER COLUMN block_size DROP DEFAULT,
   ALTER COLUMN region_features DROP DEFAULT,
   ALTER COLUMN value_type DROP DEFAULT;

ALTER TABLE ditag 
   ALTER COLUMN name DROP DEFAULT,
   ALTER COLUMN type DROP DEFAULT;

ALTER TABLE ditag_feature 
   ALTER COLUMN ditag_side DROP DEFAULT;

ALTER TABLE dna 
   ALTER COLUMN seq_region_id DROP DEFAULT;

ALTER TABLE dna_align_feature 
   ALTER COLUMN analysis_id DROP DEFAULT,
   ALTER COLUMN hit_end DROP DEFAULT,
   ALTER COLUMN hit_name DROP DEFAULT,
   ALTER COLUMN hit_start DROP DEFAULT,
   ALTER COLUMN hit_strand DROP DEFAULT,
   ALTER COLUMN seq_region_end DROP DEFAULT,
   ALTER COLUMN seq_region_id DROP DEFAULT,
   ALTER COLUMN seq_region_start DROP DEFAULT,
   ALTER COLUMN seq_region_strand DROP DEFAULT;

ALTER TABLE dnac 
   ALTER COLUMN seq_region_id DROP DEFAULT;

ALTER TABLE exon 
   ALTER COLUMN end_phase DROP DEFAULT,
   ALTER COLUMN phase DROP DEFAULT,
   ALTER COLUMN seq_region_end DROP DEFAULT,
   ALTER COLUMN seq_region_id DROP DEFAULT,
   ALTER COLUMN seq_region_start DROP DEFAULT,
   ALTER COLUMN seq_region_strand DROP DEFAULT;

ALTER TABLE exon_transcript 
   ALTER COLUMN exon_id DROP DEFAULT,
   ALTER COLUMN rank DROP DEFAULT,
   ALTER COLUMN transcript_id DROP DEFAULT;

ALTER TABLE external_synonym 
   ALTER COLUMN xref_id DROP DEFAULT;

ALTER TABLE gene 
   ALTER COLUMN analysis_id DROP DEFAULT,
   ALTER COLUMN biotype DROP DEFAULT,
   ALTER COLUMN seq_region_end DROP DEFAULT,
   ALTER COLUMN seq_region_id DROP DEFAULT,
   ALTER COLUMN seq_region_start DROP DEFAULT,
   ALTER COLUMN seq_region_strand DROP DEFAULT,
   ALTER COLUMN source DROP DEFAULT;

ALTER TABLE gene_archive 
   ALTER COLUMN gene_stable_id DROP DEFAULT,
   ALTER COLUMN mapping_session_id DROP DEFAULT,
   ALTER COLUMN transcript_stable_id DROP DEFAULT;

ALTER TABLE identity_xref 
   ALTER COLUMN object_xref_id DROP DEFAULT;

ALTER TABLE interpro 
   ALTER COLUMN id DROP DEFAULT,
   ALTER COLUMN interpro_ac DROP DEFAULT;

ALTER TABLE karyotype 
   ALTER COLUMN band DROP DEFAULT,
   ALTER COLUMN seq_region_id DROP DEFAULT,
   ALTER COLUMN stain DROP DEFAULT;

ALTER TABLE map 
   ALTER COLUMN map_name DROP DEFAULT;

ALTER TABLE mapping_session 
   ALTER COLUMN created DROP DEFAULT;

ALTER TABLE marker 
   ALTER COLUMN left_primer DROP DEFAULT,
   ALTER COLUMN max_primer_dist DROP DEFAULT,
   ALTER COLUMN min_primer_dist DROP DEFAULT,
   ALTER COLUMN right_primer DROP DEFAULT;

ALTER TABLE marker_feature 
   ALTER COLUMN analysis_id DROP DEFAULT,
   ALTER COLUMN marker_id DROP DEFAULT,
   ALTER COLUMN seq_region_end DROP DEFAULT,
   ALTER COLUMN seq_region_id DROP DEFAULT,
   ALTER COLUMN seq_region_start DROP DEFAULT;

ALTER TABLE marker_map_location 
   ALTER COLUMN chromosome_name DROP DEFAULT,
   ALTER COLUMN map_id DROP DEFAULT,
   ALTER COLUMN marker_id DROP DEFAULT,
   ALTER COLUMN marker_synonym_id DROP DEFAULT,
   ALTER COLUMN position DROP DEFAULT;

ALTER TABLE marker_synonym 
   ALTER COLUMN marker_id DROP DEFAULT;

ALTER TABLE meta 
   ALTER COLUMN meta_key DROP DEFAULT;

ALTER TABLE meta_coord 
   ALTER COLUMN coord_system_id DROP DEFAULT,
   ALTER COLUMN table_name DROP DEFAULT;

ALTER TABLE object_xref 
   ALTER COLUMN ensembl_id DROP DEFAULT,
   ALTER COLUMN xref_id DROP DEFAULT;

ALTER TABLE prediction_exon 
   ALTER COLUMN exon_rank DROP DEFAULT,
   ALTER COLUMN prediction_transcript_id DROP DEFAULT,
   ALTER COLUMN seq_region_end DROP DEFAULT,
   ALTER COLUMN seq_region_id DROP DEFAULT,
   ALTER COLUMN seq_region_start DROP DEFAULT,
   ALTER COLUMN seq_region_strand DROP DEFAULT,
   ALTER COLUMN start_phase DROP DEFAULT;

ALTER TABLE prediction_transcript 
   ALTER COLUMN analysis_id DROP DEFAULT,
   ALTER COLUMN seq_region_end DROP DEFAULT,
   ALTER COLUMN seq_region_id DROP DEFAULT,
   ALTER COLUMN seq_region_start DROP DEFAULT,
   ALTER COLUMN seq_region_strand DROP DEFAULT;

ALTER TABLE protein_align_feature 
   ALTER COLUMN analysis_id DROP DEFAULT,
   ALTER COLUMN hit_end DROP DEFAULT,
   ALTER COLUMN hit_name DROP DEFAULT,
   ALTER COLUMN hit_start DROP DEFAULT,
   ALTER COLUMN seq_region_end DROP DEFAULT,
   ALTER COLUMN seq_region_id DROP DEFAULT,
   ALTER COLUMN seq_region_start DROP DEFAULT;

ALTER TABLE protein_feature 
   ALTER COLUMN analysis_id DROP DEFAULT,
   ALTER COLUMN hit_end DROP DEFAULT,
   ALTER COLUMN hit_start DROP DEFAULT,
   ALTER COLUMN seq_end DROP DEFAULT,
   ALTER COLUMN seq_start DROP DEFAULT,
   ALTER COLUMN translation_id DROP DEFAULT;

ALTER TABLE qtl 
   ALTER COLUMN trait DROP DEFAULT;

ALTER TABLE qtl_feature 
   ALTER COLUMN analysis_id DROP DEFAULT,
   ALTER COLUMN qtl_id DROP DEFAULT,
   ALTER COLUMN seq_region_end DROP DEFAULT,
   ALTER COLUMN seq_region_id DROP DEFAULT,
   ALTER COLUMN seq_region_start DROP DEFAULT;

ALTER TABLE qtl_synonym 
   ALTER COLUMN qtl_id DROP DEFAULT,
   ALTER COLUMN source_database DROP DEFAULT,
   ALTER COLUMN source_primary_id DROP DEFAULT;

ALTER TABLE repeat_consensus 
   ALTER COLUMN repeat_class DROP DEFAULT,
   ALTER COLUMN repeat_name DROP DEFAULT,
   ALTER COLUMN repeat_type DROP DEFAULT;

ALTER TABLE repeat_feature 
   ALTER COLUMN analysis_id DROP DEFAULT,
   ALTER COLUMN repeat_consensus_id DROP DEFAULT,
   ALTER COLUMN repeat_end DROP DEFAULT,
   ALTER COLUMN repeat_start DROP DEFAULT,
   ALTER COLUMN seq_region_end DROP DEFAULT,
   ALTER COLUMN seq_region_id DROP DEFAULT,
   ALTER COLUMN seq_region_start DROP DEFAULT;

ALTER TABLE seq_region 
   ALTER COLUMN coord_system_id DROP DEFAULT,
   ALTER COLUMN name DROP DEFAULT;

ALTER TABLE simple_feature 
   ALTER COLUMN analysis_id DROP DEFAULT,
   ALTER COLUMN seq_region_end DROP DEFAULT,
   ALTER COLUMN seq_region_id DROP DEFAULT,
   ALTER COLUMN seq_region_start DROP DEFAULT,
   ALTER COLUMN seq_region_strand DROP DEFAULT;

ALTER TABLE stable_id_event 
   ALTER COLUMN type DROP DEFAULT;

ALTER TABLE transcript 
   ALTER COLUMN analysis_id DROP DEFAULT,
   ALTER COLUMN biotype DROP DEFAULT,
   ALTER COLUMN seq_region_end DROP DEFAULT,
   ALTER COLUMN seq_region_id DROP DEFAULT,
   ALTER COLUMN seq_region_start DROP DEFAULT,
   ALTER COLUMN seq_region_strand DROP DEFAULT;

ALTER TABLE translation 
   ALTER COLUMN end_exon_id DROP DEFAULT,
   ALTER COLUMN seq_end DROP DEFAULT,
   ALTER COLUMN seq_start DROP DEFAULT,
   ALTER COLUMN start_exon_id DROP DEFAULT,
   ALTER COLUMN transcript_id DROP DEFAULT;

ALTER TABLE unconventional_transcript_association 
   ALTER COLUMN gene_id DROP DEFAULT,
   ALTER COLUMN transcript_id DROP DEFAULT;

ALTER TABLE unmapped_object 
   ALTER COLUMN analysis_id DROP DEFAULT,
   ALTER COLUMN identifier DROP DEFAULT,
   ALTER COLUMN unmapped_reason_id DROP DEFAULT;

ALTER TABLE xref 
   ALTER COLUMN dbprimary_acc DROP DEFAULT,
   ALTER COLUMN display_label DROP DEFAULT;

# Patch identifier:
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_65_66_f.sql|drop_default_values');
