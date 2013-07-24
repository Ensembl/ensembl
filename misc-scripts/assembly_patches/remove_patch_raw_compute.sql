delete repeat_feature from attrib_type, seq_region_attrib, repeat_feature where attrib_type.code in ('patch_novel','patch_fix') and attrib_type.attrib_type_id = seq_region_attrib.attrib_type_id and seq_region_attrib.seq_region_id = repeat_feature.seq_region_id;

delete prediction_exon from attrib_type, seq_region_attrib, prediction_exon where attrib_type.code in ('patch_novel','patch_fix') and attrib_type.attrib_type_id = seq_region_attrib.attrib_type_id and seq_region_attrib.seq_region_id = prediction_exon.seq_region_id;

delete prediction_transcript from attrib_type, seq_region_attrib, prediction_transcript where attrib_type.code in ('patch_novel','patch_fix') and attrib_type.attrib_type_id = seq_region_attrib.attrib_type_id and seq_region_attrib.seq_region_id = prediction_transcript.seq_region_id;

delete simple_feature from attrib_type, seq_region_attrib, simple_feature where attrib_type.code in ('patch_novel','patch_fix') and attrib_type.attrib_type_id = seq_region_attrib.attrib_type_id and seq_region_attrib.seq_region_id = simple_feature.seq_region_id;

delete dna_align_feature from attrib_type, seq_region_attrib, dna_align_feature where attrib_type.code in ('patch_novel','patch_fix') and attrib_type.attrib_type_id = seq_region_attrib.attrib_type_id and seq_region_attrib.seq_region_id = dna_align_feature.seq_region_id and dna_align_feature_id not in (select feature_id from transcript_supporting_feature where feature_type = 'dna_align_feature') and dna_align_feature_id not in (select feature_id from supporting_feature where feature_type = 'dna_align_feature');

delete protein_align_feature from attrib_type, seq_region_attrib, protein_align_feature where attrib_type.code in ('patch_novel','patch_fix') and attrib_type.attrib_type_id = seq_region_attrib.attrib_type_id and seq_region_attrib.seq_region_id = protein_align_feature.seq_region_id and protein_align_feature_id not in (select feature_id from transcript_supporting_feature where feature_type = 'protein_align_feature') and protein_align_feature_id not in (select feature_id from supporting_feature where feature_type = 'protein_align_feature');

