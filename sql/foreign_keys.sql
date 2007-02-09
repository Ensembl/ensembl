# Foreign key relationships in the Ensembl schema (see table.sql)
#
# This file is intended as a reference since some of the relationships
# are not obvious.
#
# Note that these constraints are not actually used by Ensembl for 
# performance reasons, and referential integrity is enforced at the
# application level. Also MySQL currently does not support foreign
# key constraints on MyISAM tables.

ALTER table alt_allele ADD FOREIGN KEY (gene_id) REFERENCES gene(gene_id);

ALTER table analysis_description ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);

ALTER table assembly ADD FOREIGN KEY (asm_seq_region_id) REFERENCES seq_region(seq_region_id);
ALTER table assembly ADD FOREIGN KEY (cmp_seq_region_id) REFERENCES seq_region(seq_region_id);

ALTER table assembly_exception ADD FOREIGN KEY (exc_seq_region_id) REFERENCES seq_region(seq_region_id);
ALTER table assembly_exception ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);

ALTER table density_feature ADD FOREIGN KEY (density_type_id) REFERENCES density_type(density_type_id);
ALTER table density_feature ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);

ALTER table density_type ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);

ALTER table ditag_feature ADD FOREIGN KEY (ditag_id) REFERENCES ditag(ditag_id);
ALTER table ditag_feature ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);
ALTER table ditag_feature ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);

ALTER table dna ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);

ALTER table dna_align_feature ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);
ALTER table dna_align_feature ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);

ALTER table dnac ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);

ALTER table exon ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);

ALTER table exon_stable_id ADD FOREIGN KEY (exon_id) REFERENCES exon(exon_id);

ALTER table exon_transcript ADD FOREIGN KEY (exon_id) REFERENCES exon(exon_id);
ALTER table exon_transcript ADD FOREIGN KEY (transcript_id) REFERENCES transcript(transcript_id);

ALTER table external_synonym ADD FOREIGN KEY (xref_id) REFERENCES xref(xref_id);

ALTER table gene ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);
ALTER table gene ADD FOREIGN KEY (display_xref_id) REFERENCES xref(xref_id);
ALTER table gene ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);

ALTER table gene_attrib ADD FOREIGN KEY (attrib_type_id) REFERENCES attrib_type(attrib_type_id);
ALTER table gene_attrib ADD FOREIGN KEY (gene_id) REFERENCES gene(gene_id);

ALTER table gene_archive ADD FOREIGN KEY (mapping_session_id) REFERENCES mapping_session (mapping_session_id);

ALTER table gene_stable_id ADD FOREIGN KEY (gene_id) REFERENCES gene(gene_id);

ALTER table go_xref ADD FOREIGN KEY (object_xref_id) REFERENCES object_xref(object_xref_id);

ALTER table identity_xref ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);
ALTER table identity_xref ADD FOREIGN KEY (object_xref_id) REFERENCES object_xref(object_xref_id);

ALTER table karyotype ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);

ALTER table stable_id_event ADD FOREIGN KEY (mapping_session_id) REFERENCES mapping_session(mapping_session_id);

ALTER table marker_feature ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);
ALTER table marker_feature ADD FOREIGN KEY (marker_id) REFERENCES marker(marker_id);
ALTER table marker_feature ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);

ALTER table marker_map_location ADD FOREIGN KEY (map_id) REFERENCES map(map_id);
ALTER table marker_map_location ADD FOREIGN KEY (marker_id) REFERENCES marker(marker_id);
ALTER table marker_map_location ADD FOREIGN KEY (marker_synonym_id) REFERENCES marker_synonym(marker_synonym_id);

ALTER table marker_synonym ADD FOREIGN KEY (marker_id) REFERENCES marker(marker_id);

ALTER table meta_coord ADD FOREIGN KEY (coord_system_id) REFERENCES coord_system(coord_system_id);

ALTER table misc_attrib ADD FOREIGN KEY (attrib_type_id) REFERENCES attrib_type(attrib_type_id);
ALTER table misc_attrib ADD FOREIGN KEY (misc_feature_id) REFERENCES misc_feature(misc_feature_id);

ALTER table misc_feature ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);

ALTER table misc_feature_misc_set ADD FOREIGN KEY (misc_feature_id) REFERENCES misc_feature(misc_feature_id);
ALTER table misc_feature_misc_set ADD FOREIGN KEY (misc_set_id) REFERENCES misc_set(misc_set_id);

ALTER table object_xref ADD FOREIGN KEY (xref_id) REFERENCES xref(xref_id);

ALTER table oligo_array ADD FOREIGN KEY (parent_array_id) REFERENCES oligo_array(oligo_array_id);

ALTER table oligo_feature ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);
ALTER table oligo_feature ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);
ALTER table oligo_feature ADD FOREIGN KEY (oligo_probe_id) REFERENCES oligo_probe (oligo_probe_id);

ALTER table oligo_probe ADD FOREIGN KEY (oligo_array_id) REFERENCES oligo_array(oligo_array_id);

ALTER table prediction_exon ADD FOREIGN KEY (prediction_transcript_id) REFERENCES prediction_transcript(prediction_transcript_id);
ALTER table prediction_exon ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);

ALTER table prediction_transcript ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);
ALTER table prediction_transcript ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);

ALTER table protein_align_feature ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);
ALTER table protein_align_feature ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);

ALTER table protein_feature ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);
ALTER table protein_feature ADD FOREIGN KEY (translation_id) REFERENCES translation(translation_id);

ALTER table qtl ADD FOREIGN KEY (flank_marker_id_1) REFERENCES marker(marker_id);
ALTER table qtl ADD FOREIGN KEY (flank_marker_id_2) REFERENCES marker(marker_id);
ALTER table qtl ADD FOREIGN KEY (peak_marker_id) REFERENCES marker(marker_id);

ALTER table qtl_feature ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);
ALTER table qtl_feature ADD FOREIGN KEY (qtl_id) REFERENCES qtl(qtl_id);
ALTER table qtl_feature ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);

ALTER table qtl_synonym ADD FOREIGN KEY (qtl_id) REFERENCES qtl(qtl_id);

ALTER table regulatory_factor_coding ADD FOREIGN KEY (gene_id) REFERENCES gene(gene_id);
ALTER table regulatory_factor_coding ADD FOREIGN KEY (regulatory_factor_id) REFERENCES regulatory_factor(regulatory_factor_id);
ALTER table regulatory_factor_coding ADD FOREIGN KEY (transcript_id) REFERENCES transcript(transcript_id);

ALTER table regulatory_feature ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);
ALTER table regulatory_feature ADD FOREIGN KEY (regulatory_factor_id) REFERENCES regulatory_factor(regulatory_factor_id);
ALTER table regulatory_feature ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);

ALTER table regulatory_search_region ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);
ALTER table regulatory_search_region ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);

ALTER table repeat_feature ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);
ALTER table repeat_feature ADD FOREIGN KEY (repeat_consensus_id) REFERENCES repeat_consensus(repeat_consensus_id);
ALTER table repeat_feature ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);

ALTER table seq_region ADD FOREIGN KEY (coord_system_id) REFERENCES coord_system(coord_system_id);

ALTER table seq_region_attrib ADD FOREIGN KEY (attrib_type_id) REFERENCES attrib_type(attrib_type_id);
ALTER table seq_region_attrib ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);

ALTER table simple_feature ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);
ALTER table simple_feature ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);

ALTER table supporting_feature ADD FOREIGN KEY (exon_id) REFERENCES exon(exon_id);

ALTER table transcript ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);
ALTER table transcript ADD FOREIGN KEY (display_xref_id) REFERENCES xref(xref_id);
ALTER table transcript ADD FOREIGN KEY (gene_id) REFERENCES gene(gene_id);
ALTER table transcript ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);

ALTER table transcript_attrib ADD FOREIGN KEY (attrib_type_id) REFERENCES attrib_type(attrib_type_id);
ALTER table transcript_attrib ADD FOREIGN KEY (transcript_id) REFERENCES transcript(transcript_id);

ALTER table transcript_stable_id ADD FOREIGN KEY (transcript_id) REFERENCES transcript(transcript_id);

ALTER table transcript_supporting_feature ADD FOREIGN KEY (transcript_id) REFERENCES transcript(transcript_id);

ALTER table translation ADD FOREIGN KEY (end_exon_id) REFERENCES exon(exon_id);
ALTER table translation ADD FOREIGN KEY (start_exon_id) REFERENCES exon(exon_id);
ALTER table translation ADD FOREIGN KEY (transcript_id) REFERENCES transcript(transcript_id);

ALTER table translation_attrib ADD FOREIGN KEY (attrib_type_id) REFERENCES attrib_type(attrib_type_id);
ALTER table translation_attrib ADD FOREIGN KEY (translation_id) REFERENCES translation(translation_id);

ALTER table translation_stable_id ADD FOREIGN KEY (translation_id) REFERENCES translation(translation_id);

ALTER table unconventional_transcript_association ADD FOREIGN KEY (gene_id) REFERENCES gene(gene_id);
ALTER table unconventional_transcript_association ADD FOREIGN KEY (transcript_id) REFERENCES transcript(transcript_id);

ALTER table unmapped_object ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);
#ALTER table unmapped_object ADD FOREIGN KEY (external_db_id) REFERENCES external_db(external_db_id);
ALTER table unmapped_object ADD FOREIGN KEY (unmapped_reason_id) REFERENCES unmapped_reason(unmapped_reason_id);

ALTER table xref ADD FOREIGN KEY (external_db_id) REFERENCES external_db(external_db_id);

