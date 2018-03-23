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

# Foreign key relationships in the Ensembl schema (see table.sql)
#
# This file is intended as a reference since some of the relationships
# are not obvious.
#
# Note that these constraints are not actually used by Ensembl for 
# performance reasons, and referential integrity is enforced at the
# application level. Also MySQL currently does not support foreign
# key constraints on MyISAM tables.
#
# Foreign keys are grouped by category

# Assembly

ALTER TABLE karyotype ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region (seq_region_id);

ALTER TABLE assembly ADD FOREIGN KEY (asm_seq_region_id) REFERENCES seq_region (seq_region_id);
ALTER TABLE assembly ADD FOREIGN KEY (cmp_seq_region_id) REFERENCES seq_region (seq_region_id);

ALTER TABLE assembly_exception ADD FOREIGN KEY (exc_seq_region_id) REFERENCES seq_region (seq_region_id);
ALTER TABLE assembly_exception ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region (seq_region_id);

ALTER TABLE data_file ADD FOREIGN KEY (coord_system_id) REFERENCES coord_system (coord_system_id);
ALTER TABLE data_file ADD FOREIGN KEY (analysis_id) REFERENCES analysis (analysis_id);

ALTER TABLE dna ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region (seq_region_id);

ALTER TABLE meta ADD FOREIGN KEY (species_id) REFERENCES coord_system (species_id);

ALTER TABLE meta_coord ADD FOREIGN KEY (coord_system_id) REFERENCES coord_system (coord_system_id);

ALTER TABLE seq_region ADD FOREIGN KEY (coord_system_id) REFERENCES coord_system (coord_system_id);

ALTER TABLE seq_region_synonym ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region (seq_region_id);

ALTER TABLE seq_region_attrib ADD FOREIGN KEY (attrib_type_id) REFERENCES attrib_type (attrib_type_id);
ALTER TABLE seq_region_attrib ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region (seq_region_id);

# Fundamental

ALTER TABLE alt_allele ADD FOREIGN KEY (gene_id) REFERENCES gene (gene_id);
ALTER TABLE alt_allele ADD FOREIGN KEY (alt_allele_group_id) REFERENCES alt_allele_group (alt_allele_group_id);
ALTER TABLE alt_allele_attrib ADD FOREIGN KEY (alt_allele_id) REFERENCES alt_allele (alt_allele_id);

ALTER TABLE analysis_description ADD FOREIGN KEY (analysis_id) REFERENCES analysis (analysis_id);

ALTER TABLE dna_align_feature ADD FOREIGN KEY (analysis_id) REFERENCES analysis (analysis_id);
ALTER TABLE dna_align_feature ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region (seq_region_id);
ALTER TABLE dna_align_feature ADD FOREIGN KEY (external_db_id) REFERENCES external_db (external_db_id);

ALTER TABLE exon ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region (seq_region_id);

ALTER TABLE exon_transcript ADD FOREIGN KEY (exon_id) REFERENCES exon (exon_id);
ALTER TABLE exon_transcript ADD FOREIGN KEY (transcript_id) REFERENCES transcript (transcript_id);

ALTER TABLE gene ADD FOREIGN KEY (analysis_id) REFERENCES analysis (analysis_id);
ALTER TABLE gene ADD FOREIGN KEY (display_xref_id) REFERENCES xref (xref_id);
ALTER TABLE gene ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region (seq_region_id);
ALTER TABLE gene ADD FOREIGN KEY (canonical_transcript_id) REFERENCES transcript (transcript_id);

ALTER TABLE gene_attrib ADD FOREIGN KEY (attrib_type_id) REFERENCES attrib_type (attrib_type_id);
ALTER TABLE gene_attrib ADD FOREIGN KEY (gene_id) REFERENCES gene (gene_id);

ALTER TABLE protein_align_feature ADD FOREIGN KEY (analysis_id) REFERENCES analysis (analysis_id);
ALTER TABLE protein_align_feature ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region (seq_region_id);
ALTER TABLE protein_align_feature ADD FOREIGN KEY (external_db_id) REFERENCES external_db (external_db_id);

ALTER TABLE protein_feature ADD FOREIGN KEY (analysis_id) REFERENCES analysis (analysis_id);
ALTER TABLE protein_feature ADD FOREIGN KEY (translation_id) REFERENCES translation (translation_id);

ALTER TABLE supporting_feature ADD FOREIGN KEY (exon_id) REFERENCES exon (exon_id);

ALTER TABLE transcript ADD FOREIGN KEY (analysis_id) REFERENCES analysis (analysis_id);
ALTER TABLE transcript ADD FOREIGN KEY (display_xref_id) REFERENCES xref (xref_id);
ALTER TABLE transcript ADD FOREIGN KEY (gene_id) REFERENCES gene (gene_id);
ALTER TABLE transcript ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region (seq_region_id);
ALTER TABLE transcript ADD FOREIGN KEY (canonical_translation_id) REFERENCES translation (translation_id);

ALTER TABLE transcript_attrib ADD FOREIGN KEY (attrib_type_id) REFERENCES attrib_type (attrib_type_id);
ALTER TABLE transcript_attrib ADD FOREIGN KEY (transcript_id) REFERENCES transcript (transcript_id);

ALTER TABLE transcript_supporting_feature ADD FOREIGN KEY (transcript_id) REFERENCES transcript (transcript_id);

ALTER TABLE translation ADD FOREIGN KEY (end_exon_id) REFERENCES exon (exon_id);
ALTER TABLE translation ADD FOREIGN KEY (start_exon_id) REFERENCES exon (exon_id);
ALTER TABLE translation ADD FOREIGN KEY (transcript_id) REFERENCES transcript (transcript_id);

ALTER TABLE translation_attrib ADD FOREIGN KEY (attrib_type_id) REFERENCES attrib_type (attrib_type_id);
ALTER TABLE translation_attrib ADD FOREIGN KEY (translation_id) REFERENCES translation (translation_id);

# Analysis

ALTER TABLE density_feature ADD FOREIGN KEY (density_type_id) REFERENCES density_type (density_type_id);
ALTER TABLE density_feature ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region (seq_region_id);

ALTER TABLE density_type ADD FOREIGN KEY (analysis_id) REFERENCES analysis (analysis_id);

ALTER TABLE ditag_feature ADD FOREIGN KEY (ditag_id) REFERENCES ditag (ditag_id);
ALTER TABLE ditag_feature ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region (seq_region_id);
ALTER TABLE ditag_feature ADD FOREIGN KEY (analysis_id) REFERENCES analysis (analysis_id);

ALTER TABLE intron_supporting_evidence ADD FOREIGN KEY (analysis_id) REFERENCES analysis (analysis_id);
ALTER TABLE intron_supporting_evidence ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region (seq_region_id);

ALTER TABLE marker ADD FOREIGN KEY (display_marker_synonym_id) REFERENCES marker_synonym  (marker_synonym_id);

ALTER TABLE marker_feature ADD FOREIGN KEY (analysis_id) REFERENCES analysis (analysis_id);
ALTER TABLE marker_feature ADD FOREIGN KEY (marker_id) REFERENCES marker (marker_id);
ALTER TABLE marker_feature ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region (seq_region_id);

ALTER TABLE marker_map_location ADD FOREIGN KEY (map_id) REFERENCES map (map_id);
ALTER TABLE marker_map_location ADD FOREIGN KEY (marker_id) REFERENCES marker (marker_id);
ALTER TABLE marker_map_location ADD FOREIGN KEY (marker_synonym_id) REFERENCES marker_synonym (marker_synonym_id);

ALTER TABLE marker_synonym ADD FOREIGN KEY (marker_id) REFERENCES marker (marker_id);

ALTER TABLE misc_attrib ADD FOREIGN KEY (attrib_type_id) REFERENCES attrib_type (attrib_type_id);
ALTER TABLE misc_attrib ADD FOREIGN KEY (misc_feature_id) REFERENCES misc_feature (misc_feature_id);

ALTER TABLE misc_feature ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region (seq_region_id);

ALTER TABLE misc_feature_misc_set ADD FOREIGN KEY (misc_feature_id) REFERENCES misc_feature (misc_feature_id);
ALTER TABLE misc_feature_misc_set ADD FOREIGN KEY (misc_set_id) REFERENCES misc_set (misc_set_id);

ALTER TABLE simple_feature ADD FOREIGN KEY (analysis_id) REFERENCES analysis (analysis_id);
ALTER TABLE simple_feature ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region (seq_region_id);

ALTER TABLE prediction_exon ADD FOREIGN KEY (prediction_transcript_id) REFERENCES prediction_transcript (prediction_transcript_id);
ALTER TABLE prediction_exon ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region (seq_region_id);

ALTER TABLE prediction_transcript ADD FOREIGN KEY (analysis_id) REFERENCES analysis (analysis_id);
ALTER TABLE prediction_transcript ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region (seq_region_id);

ALTER TABLE repeat_feature ADD FOREIGN KEY (analysis_id) REFERENCES analysis (analysis_id);
ALTER TABLE repeat_feature ADD FOREIGN KEY (repeat_consensus_id) REFERENCES repeat_consensus (repeat_consensus_id);
ALTER TABLE repeat_feature ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region (seq_region_id);

# Xrefs

ALTER TABLE associated_xref ADD FOREIGN KEY (associated_group_id) REFERENCES associated_group (associated_group_id);
ALTER TABLE associated_xref ADD FOREIGN KEY (object_xref_id) REFERENCES object_xref (object_xref_id);
ALTER TABLE associated_xref ADD FOREIGN KEY (xref_id) REFERENCES xref (xref_id);

ALTER TABLE dependent_xref ADD FOREIGN KEY (object_xref_id) REFERENCES object_xref (object_xref_id);
ALTER TABLE dependent_xref ADD FOREIGN KEY (master_xref_id) REFERENCES xref (xref_id);
ALTER TABLE dependent_xref ADD FOREIGN KEY (dependent_xref_id) REFERENCES xref (xref_id);

ALTER TABLE external_synonym ADD FOREIGN KEY (xref_id) REFERENCES xref (xref_id);

ALTER TABLE identity_xref ADD FOREIGN KEY (object_xref_id) REFERENCES object_xref (object_xref_id);

ALTER TABLE object_xref ADD FOREIGN KEY (xref_id) REFERENCES xref (xref_id);
ALTER TABLE object_xref ADD FOREIGN KEY (analysis_id) REFERENCES analysis (analysis_id);

ALTER TABLE ontology_xref ADD FOREIGN KEY (object_xref_id) REFERENCES object_xref (object_xref_id);
ALTER TABLE ontology_xref ADD FOREIGN KEY (source_xref_id) REFERENCES xref (xref_id);

ALTER TABLE unmapped_object ADD FOREIGN KEY (analysis_id) REFERENCES analysis (analysis_id);
ALTER TABLE unmapped_object ADD FOREIGN KEY (external_db_id) REFERENCES external_db (external_db_id);
ALTER TABLE unmapped_object ADD FOREIGN KEY (unmapped_reason_id) REFERENCES unmapped_reason (unmapped_reason_id);

ALTER TABLE xref ADD FOREIGN KEY (external_db_id) REFERENCES external_db (external_db_id);

# Mapping

ALTER TABLE gene_archive ADD FOREIGN KEY (mapping_session_id) REFERENCES mapping_session (mapping_session_id);
ALTER TABLE gene_archive ADD FOREIGN KEY (peptide_archive_id) REFERENCES peptide_archive (peptide_archive_id);

ALTER TABLE stable_id_event ADD FOREIGN KEY (mapping_session_id) REFERENCES mapping_session (mapping_session_id);

ALTER TABLE seq_region_mapping ADD FOREIGN KEY (internal_seq_region_id) REFERENCES seq_region (seq_region_id);
ALTER TABLE seq_region_mapping ADD FOREIGN KEY (mapping_set_id) REFERENCES mapping_set (mapping_set_id);

# Misc

ALTER TABLE operon ADD FOREIGN KEY (analysis_id) REFERENCES analysis (analysis_id);
ALTER TABLE operon ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region (seq_region_id);

ALTER TABLE operon_transcript ADD FOREIGN KEY (analysis_id) REFERENCES analysis (analysis_id);
ALTER TABLE operon_transcript ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region (seq_region_id); 
ALTER TABLE operon_transcript ADD FOREIGN KEY (operon_id) REFERENCES operon (operon_id);

ALTER TABLE operon_transcript_gene ADD FOREIGN KEY (operon_transcript_id) REFERENCES operon_transcript (operon_transcript_id);
ALTER TABLE operon_transcript_gene ADD FOREIGN KEY (gene_id) REFERENCES gene (gene_id);









