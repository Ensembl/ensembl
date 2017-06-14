Ensembl Deprecated Methods
===================

This file contains the list of methods deprecated in the Ensembl core API.
A method is deprecated when it is not functional any more (schema/data change) or has been replaced by a better one.
Backwards compatibility is provided whenever possible.
When a method is deprecated, a deprecation warning is thrown whenever the method is used.
The warning also contains instructions on replacing the deprecated method and when it will be removed.
A year after deprecation (4 Ensembl releases), the method is removed from the API.

### Removed in Ensembl Release 91 ###

 - Bio::EnsEMBL::DBSQL::**BaseAdaptor**::*dump_data()*
 - Bio::EnsEMBL::DBSQL::**BaseAdaptor**::*get_dumped_data()*

### Removed in Ensembl Release 90 ###

 - Bio::EnsEMBL::**Gene**::*is_known()*
 - Bio::EnsEMBL::**Gene**::*status()*
 - Bio::EnsEMBL::**Transcript**::*is_known()*
 - Bio::EnsEMBL::**Transcript**::*status()*

### Removed in Ensembl Release 88 ###

 - Bio::EnsEMBL::**Slice**::*get_all_VariationFeatures()*
 - Bio::EnsEMBL::**Slice**::*get_all_somatic_VariationFeatures()*
 - Bio::EnsEMBL::**Slice**::*get_all_somatic_VariationFeatures_by_source()*
 - Bio::EnsEMBL::**Slice**::*get_all_VariationFeatures_with_phenotype()*
 - Bio::EnsEMBL::**Slice**::*get_all_somatic_VariationFeatures_with_phenotype()*
 - Bio::EnsEMBL::**Slice**::*get_all_StructuralVariationFeatures()*
 - Bio::EnsEMBL::**Slice**::*get_all_somatic_StructuralVariationFeatures()*
 - Bio::EnsEMBL::**Slice**::*get_all_StructuralVariationFeatures_by_size_range()*
 - Bio::EnsEMBL::**Slice**::*get_all_somatic_StructuralVariationFeatures_by_size_range()*
 - Bio::EnsEMBL::**Slice**::*get_all_StructuralVariationFeatures_by_Study()*
 - Bio::EnsEMBL::**Slice**::*get_all_StructuralVariationFeatures_by_source()*
 - Bio::EnsEMBL::**Slice**::*get_all_somatic_StructuralVariationFeatures_by_source()*
 - Bio::EnsEMBL::**Slice**::*get_all_VariationFeatures_by_VariationSet()*
 - Bio::EnsEMBL::**Slice**::*get_all_StructuralVariationFeatures_by_VariationSet()*
 - Bio::EnsEMBL::**Slice**::*get_all_PhenotypeFeatures()*
 - Bio::EnsEMBL::**Slice**::*get_all_CopyNumberVariantProbeFeatures()*
 - Bio::EnsEMBL::**Slice**::*get_all_genotyped_VariationFeatures()*
 - Bio::EnsEMBL::**Slice**::*get_all_VariationFeatures_by_Population()*
 - Bio::EnsEMBL::**Slice**::*get_by_strain()*
 - Bio::EnsEMBL::**StrainSlice**::*remove_indels()*
 - Bio::EnsEMBL::**StrainSlice**::*get_original_seq_region_position()*
 - Bio::EnsEMBL::**StrainSlice**::*get_all_VariationFeatures()*
 - Bio::EnsEMBL::**StrainSlice**::*mapper()*
 - Bio::EnsEMBL::**StrainSlice**::*subseq()*
 - Bio::EnsEMBL::**StrainSlice**::*ref_subseq()*
 - Bio::EnsEMBL::**StrainSlice**::*sub_Slice()*
 - Bio::EnsEMBL::**StrainSlice**::*_convert_difference()*
 - Bio::EnsEMBL::**StrainSlice**::*get_all_differences_StrainSlice()*
 - Bio::EnsEMBL::**StrainSlice**::*get_all_AlleleFeatures_Slice()*
 - Bio::EnsEMBL::**StrainSlice**::*get_all_AlleleFeature()*
 - Bio::EnsEMBL::**StrainSlice**::*_add_coverage_information()*
 - Bio::EnsEMBL::**StrainSlice**::*expanded_length()*
 - Bio::EnsEMBL::**StrainSlice**::*seq()*
 - Bio::EnsEMBL::**StrainSlice**::*display_Slice_name()*
 - Bio::EnsEMBL::**StrainSlice**::*sample()*
 - Bio::EnsEMBL::**StrainSlice**::*strain_name()*
 - Bio::EnsEMBL::**StrainSlice**::*_filter_af_by_coverage()*
 - Bio::EnsEMBL::**StrainSlice**::*new()*
 - Bio::EnsEMBL::DBSQL::**StrainSliceAdaptor**::*new()*
 - Bio::EnsEMBL::DBSQL::**StrainSliceAdaptor**::*fetch_by_name()*

### Removed in Ensembl Release 87 ###

 - Bio::EnsEMBL::**AssemblyMapper**::*in_assembly()*
 - Bio::EnsEMBL::**AssemblyMapper**::*map_coordinates_to_assembly()*
 - Bio::EnsEMBL::**AssemblyMapper**::*fast_to_assembly()*
 - Bio::EnsEMBL::**AssemblyMapper**::*map_coordinates_to_rawcontig()*
 - Bio::EnsEMBL::**AssemblyMapper**::*list_contig_ids()*
 - Bio::EnsEMBL::**ChainedAssemblyMapper**::*in_assembly()*
 - Bio::EnsEMBL::**ChainedAssemblyMapper**::*map_coordinates_to_assembly()*
 - Bio::EnsEMBL::**ChainedAssemblyMapper**::*fast_to_assembly()*
 - Bio::EnsEMBL::**ChainedAssemblyMapper**::*map_coordinates_to_rawcontig()*
 - Bio::EnsEMBL::**ChainedAssemblyMapper**::*list_contig_ids()*
 - Bio::EnsEMBL::**DBEntry**::*get_synonyms()*
 - Bio::EnsEMBL::DBSQL::**AssemblyMapperAdaptor**::*register_region()*
 - Bio::EnsEMBL::DBSQL::**AssemblyMapperAdaptor**::*register_contig()*
 - Bio::EnsEMBL::DBSQL::**AssemblyMapperAdaptor**::*fetch_by_type()*
 - Bio::EnsEMBL::DBSQL::**KaryotypeBandAdaptor**::*fetch_by_chr_band()*
 - Bio::EnsEMBL::DBSQL::**TranslationAdaptor**::*fetch_all_by_**DBEntry**()*
 - Bio::EnsEMBL::DBSQL::**TranslationAdaptor**::*get_stable_entry_info()*
 - Bio::EnsEMBL::DBSQL::**AltAlleleGroupAdaptor**::*fetch_all_Groups()*
 - Bio::EnsEMBL::DBSQL::**AltAlleleGroupAdaptor**::*fetch_all_Groups_by_type()*
 - Bio::EnsEMBL::DBSQL::**AltAlleleGroupAdaptor**::*fetch_Group_by_id()*
 - Bio::EnsEMBL::DBSQL::**AltAlleleGroupAdaptor**::*fetch_Group_by_Gene_dbID()*
 - Bio::EnsEMBL::DBSQL::**AnalysisAdaptor**::*feature_classes()*
 - Bio::EnsEMBL::DBSQL::**BaseAlignFeatureAdaptor**::*fetch_all_by_RawContig_and_pid()*
 - Bio::EnsEMBL::DBSQL::**BaseFeatureAdaptor**::*fetch_all_by_RawContig_constraint()*
 - Bio::EnsEMBL::DBSQL::**BaseFeatureAdaptor**::*fetch_all_by_RawContig()*
 - Bio::EnsEMBL::DBSQL::**BaseFeatureAdaptor**::*fetch_all_by_RawContig_and_score()*
 - Bio::EnsEMBL::DBSQL::**BaseFeatureAdaptor**::*remove_by_RawContig()*
 - Bio::EnsEMBL::DBSQL::**DBAdaptor**::*db_handle()*
 - Bio::EnsEMBL::DBSQL::**DBAdaptor**::*port()*
 - Bio::EnsEMBL::DBSQL::**DBAdaptor**::*driver()*
 - Bio::EnsEMBL::DBSQL::**DBAdaptor**::*password()*
 - Bio::EnsEMBL::DBSQL::**DBAdaptor**::*username()*
 - Bio::EnsEMBL::DBSQL::**DBAdaptor**::*host()*
 - Bio::EnsEMBL::DBSQL::**DBAdaptor**::*reconnect_when_lost()*
 - Bio::EnsEMBL::DBSQL::**DBAdaptor**::*disconnect_when_inactive()*
 - Bio::EnsEMBL::DBSQL::**DBAdaptor**::*dbname()*
 - Bio::EnsEMBL::DBSQL::**DBAdaptor**::*prepare()*
 - Bio::EnsEMBL::DBSQL::**DBAdaptor**::*list_supported_assemblies()*
 - Bio::EnsEMBL::DBSQL::**DBAdaptor**::*assembly_type()*
 - Bio::EnsEMBL::DBSQL::**DBAdaptor**::*db()*
 - Bio::EnsEMBL::DBSQL::**DBConnection**::*group()*
 - Bio::EnsEMBL::DBSQL::**DBConnection**::*species()*
 - Bio::EnsEMBL::DBSQL::**DBEntryAdaptor**::*geneids_by_extids()*
 - Bio::EnsEMBL::DBSQL::**DBEntryAdaptor**::*translationids_by_extids()*
 - Bio::EnsEMBL::DBSQL::**DBEntryAdaptor**::*transcriptids_by_extids()*
 - Bio::EnsEMBL::DBSQL::**DataFileAdaptor**::*DataFile_to_extension()*
 - Bio::EnsEMBL::DBSQL::**ExonAdaptor**::*get_stable_entry_info()*
 - Bio::EnsEMBL::DBSQL::**ExonAdaptor**::*fetch_all_by_gene_id()*
 - Bio::EnsEMBL::DBSQL::**GeneAdaptor**::*fetch_nearest_Gene_by_Feature()*
 - Bio::EnsEMBL::DBSQL::**GeneAdaptor**::*fetch_by_maximum_DBLink()*
 - Bio::EnsEMBL::DBSQL::**GeneAdaptor**::*get_display_xref()*
 - Bio::EnsEMBL::DBSQL::**GeneAdaptor**::*get_description()*
 - Bio::EnsEMBL::DBSQL::**GeneAdaptor**::*fetch_all_by_**DBEntry**()*
 - Bio::EnsEMBL::DBSQL::**GeneAdaptor**::*get_stable_entry_info()*
 - Bio::EnsEMBL::DBSQL::**GeneAdaptor**::*fetch_by_Peptide_id()*
 - Bio::EnsEMBL::DBSQL::**MetaContainer**::*get_Species()*
 - Bio::EnsEMBL::DBSQL::**MetaContainer**::*get_default_assembly()*
 - Bio::EnsEMBL::DBSQL::**ProteinFeatureAdaptor**::*fetch_by_translation_id()*
 - Bio::EnsEMBL::DBSQL::**ProteinFeatureAdaptor**::*fetch_all_by_feature_and_dbID()*
 - Bio::EnsEMBL::DBSQL::**RepeatConsensusAdaptor**::*fetch_by_class_seq()*
 - Bio::EnsEMBL::DBSQL::**SliceAdaptor**::*fetch_by_mapfrag()*
 - Bio::EnsEMBL::DBSQL::**SliceAdaptor**::*fetch_by_chr_start_end()*
 - Bio::EnsEMBL::DBSQL::**SliceAdaptor**::*fetch_by_contig_name()*
 - Bio::EnsEMBL::DBSQL::**SliceAdaptor**::*fetch_by_clone_accession()*
 - Bio::EnsEMBL::DBSQL::**SliceAdaptor**::*fetch_by_supercontig_name()*
 - Bio::EnsEMBL::DBSQL::**SliceAdaptor**::*list_overlapping_supercontigs()*
 - Bio::EnsEMBL::DBSQL::**SliceAdaptor**::*fetch_by_chr_name()*
 - Bio::EnsEMBL::DBSQL::**TranscriptAdaptor**::*get_display_xref()*
 - Bio::EnsEMBL::DBSQL::**TranscriptAdaptor**::*get_stable_entry_info()*
 - Bio::EnsEMBL::DBSQL::**TranscriptAdaptor**::*fetch_all_by_DBEntry()*
 - Bio::EnsEMBL::DBSQL::**SequenceAdaptor**::*fetch_by_assembly_location()*
 - Bio::EnsEMBL::DBSQL::**SequenceAdaptor**::*fetch_by_RawContig_start_end_strand()*
 - Bio::EnsEMBL::**Exon**::*temporary_id()*
 - Bio::EnsEMBL::**Exon**::*created()*
 - Bio::EnsEMBL::**Exon**::*modified()*
 - Bio::EnsEMBL::**Exon**::*type()*
 - Bio::EnsEMBL::**FeaturePair**::*feature1()*
 - Bio::EnsEMBL::**FeaturePair**::*feature2()*
 - Bio::EnsEMBL::**FeaturePair**::*set_featurepair_fields()*
 - Bio::EnsEMBL::**FeaturePair**::*gffstring()*
 - Bio::EnsEMBL::**FeaturePair**::*hphase()*
 - Bio::EnsEMBL::**FeaturePair**::*hend_phase()*
 - Bio::EnsEMBL::**Feature**::*contig()*
 - Bio::EnsEMBL::**Feature**::*id()*
 - Bio::EnsEMBL::**Gene**::*add_DBLink()*
 - Bio::EnsEMBL::**Gene**::*temporary_id()*
 - Bio::EnsEMBL::**Gene**::*chr_name()*
 - Bio::EnsEMBL::**Gene**::*type()*
 - Bio::EnsEMBL::**Gene**::*confidence()*
 - Bio::EnsEMBL::**IdentityXref**::*query_identity()*
 - Bio::EnsEMBL::**IdentityXref**::*target_identity()*
 - Bio::EnsEMBL::**IdentityXref**::*translation_start()*
 - Bio::EnsEMBL::**IdentityXref**::*translation_end()*
 - Bio::EnsEMBL::**IdentityXref**::*query_start()*
 - Bio::EnsEMBL::**IdentityXref**::*query_end()*
 - Bio::EnsEMBL::**KaryotypeBand**::*chr_name()*
 - Bio::EnsEMBL::Map::DBSQL::**MarkerFeatureAdaptor**::*fetch_all_by_RawContig_and_priority()*
 - Bio::EnsEMBL::Map::**DitagFeature**::*fetch_ditag()*
 - Bio::EnsEMBL::Map::**MapLocation**::*chromosome()*
 - Bio::EnsEMBL::**OperonTranscript**::*add_gene()*
 - Bio::EnsEMBL::**PredictionTranscript**::*get_exon_count()*
 - Bio::EnsEMBL::**PredictionTranscript**::*get_cdna()*
 - Bio::EnsEMBL::**Registry**::*load_registry_with_web_adaptors()*
 - Bio::EnsEMBL::**Root**::*throw()*
 - Bio::EnsEMBL::**Root**::*warn()*
 - Bio::EnsEMBL::**Root**::*verbose()*
 - Bio::EnsEMBL::**Root**::*stack_trace_dump()*
 - Bio::EnsEMBL::**Root**::*stack_trace()*
 - Bio::EnsEMBL::**Slice**::*get_all_SNPs()*
 - Bio::EnsEMBL::**Slice**::*get_all_genotyped_SNPs()*
 - Bio::EnsEMBL::**Slice**::*get_all_supercontig_Slices()*
 - Bio::EnsEMBL::**Slice**::*get_Chromosome()*
 - Bio::EnsEMBL::**Slice**::*chr_name()*
 - Bio::EnsEMBL::**Slice**::*chr_start()*
 - Bio::EnsEMBL::**Slice**::*chr_end()*
 - Bio::EnsEMBL::**Slice**::*assembly_type()*
 - Bio::EnsEMBL::**Slice**::*dbID()*
 - Bio::EnsEMBL::**Slice**::*get_all_MapFrags()*
 - Bio::EnsEMBL::**Slice**::*has_MapSet()*
 - Bio::EnsEMBL::**StrainSlice**::*get_all_differences_Slice()*
 - Bio::EnsEMBL::**Transcript**::*created()*
 - Bio::EnsEMBL::**Transcript**::*modified()*
 - Bio::EnsEMBL::**Transcript**::*temporary_id()*
 - Bio::EnsEMBL::**Transcript**::*type()*
 - Bio::EnsEMBL::**Transcript**::*confidence()*
 - Bio::EnsEMBL::**Translation**::*temporary_id()*
 - Bio::EnsEMBL::Utils::**ConversionSupport**::*user_confirm()*

 
### Removed in Ensembl Release 84 ###
 - Bio::EnsEMBL::DBSQL::**CoordSystemAdaptor**::*_fetch_by_attrib()*
 - Bio::EnsEMBL::DBSQL::**CoordSystemAdaptor**::*_fetch_all_by_attrib()*
 - Bio::EnsEMBL::DBSQL::**DBAdaptor**::*source()*
 - Bio::EnsEMBL::DBSQL::**SliceAdaptor**::*fetch_by_band()*
 - Bio::EnsEMBL::DBSQL::**MetaContainer**::*get_short_name()*
 - Bio::EnsEMBL::DBSQL::**MetaContainer**::*get_max_assembly_contig()*
 - Bio::EnsEMBL::**DBEntry**::*ensembl_object_type()*
 - Bio::EnsEMBL::**DBEntry**::*ensembl_id()*
 - Bio::EnsEMBL::**Exon**::*_get_stable_entry_info()*
 - Bio::EnsEMBL::**FeaturePair**::*validate()*
 - Bio::EnsEMBL::**FeaturePair**::*validate_prot_feature()*
 - Bio::EnsEMBL::**PredictionTranscript**::*set_exon_count()*
 - Bio::EnsEMBL::**Root**::*_rearrange()*
 - Bio::EnsEMBL::**SeqFeatureI**::*analysis()*
 - Bio::EnsEMBL::**SeqFeatureI**::*validate()*
 - Bio::EnsEMBL::**SeqFeatureI**::*id()*
 - Bio::EnsEMBL::**SeqFeatureI**::*percent_id()*
 - Bio::EnsEMBL::**SeqFeatureI**::*e_value()*
 - Bio::EnsEMBL::**SeqFeatureI**::*phase()*
 - Bio::EnsEMBL::**SeqFeatureI**::*end_phase()*
 - Bio::EnsEMBL::**SeqFeatureI**::*location()*
 - Bio::EnsEMBL::**Slice**::*get_all_SNPs_transcripts()*
 - Bio::EnsEMBL::**Slice**::*get_all_AffyFeatures()*
 - Bio::EnsEMBL::**Slice**::*get_all_OligoFeatures()*
 - Bio::EnsEMBL::**Slice**::*get_all_OligoFeatures_by_type()*
 - Bio::EnsEMBL::**Slice**::*get_tiling_path()*
 - Bio::EnsEMBL::**Transcript**::*sort()*
 - Bio::EnsEMBL::**Transcript**::*_translation_id()*
