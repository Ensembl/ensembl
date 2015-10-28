This file contains the list of methods deprecated in the Ensembl core API.
A method is deprecated when it is not functional any more (schema/data change) or has been replaced by a better one.
Backwards compatibility is provided whenever possible.
When a method is deprecated, a deprecation warning is thrown whenever the method is used.
The warning also contains instructions on replacing the deprecated method and when it will be removed.
A year after deprecation (4 Ensembl releases), the method is removed from the API.

Removed in e87
Bio::EnsEMBL::AssemblyMapper::in_assembly
Bio::EnsEMBL::AssemblyMapper::map_coordinates_to_assembly
Bio::EnsEMBL::AssemblyMapper::fast_to_assembly
Bio::EnsEMBL::AssemblyMapper::map_coordinates_to_rawcontig
Bio::EnsEMBL::AssemblyMapper::list_contig_ids
Bio::EnsEMBL::ChainedAssemblyMapper::in_assembly
Bio::EnsEMBL::ChainedAssemblyMapper::map_coordinates_to_assembly
Bio::EnsEMBL::ChainedAssemblyMapper::fast_to_assembly
Bio::EnsEMBL::ChainedAssemblyMapper::map_coordinates_to_rawcontig
Bio::EnsEMBL::ChainedAssemblyMapper::list_contig_ids
Bio::EnsEMBL::DBEntry::get_synonyms
Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor::register_region
Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor::register_contig
Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor::fetch_by_type
Bio::EnsEMBL::DBSQL::KaryotypeBandAdaptor::fetch_by_chr_band
Bio::EnsEMBL::DBSQL::TranslationAdaptor::fetch_all_by_DBEntry
Bio::EnsEMBL::DBSQL::TranslationAdaptor::get_stable_entry_info
Bio::EnsEMBL::DBSQL::AltAlleleGroupAdaptor::fetch_all_Groups
Bio::EnsEMBL::DBSQL::AltAlleleGroupAdaptor::fetch_all_Groups_by_type
Bio::EnsEMBL::DBSQL::AltAlleleGroupAdaptor::fetch_Group_by_id
Bio::EnsEMBL::DBSQL::AltAlleleGroupAdaptor::fetch_Group_by_Gene_dbID
Bio::EnsEMBL::DBSQL::AnalysisAdaptor::feature_classes
Bio::EnsEMBL::DBSQL::BaseAlignFeatureAdaptor::fetch_all_by_RawContig_and_pid
Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor::fetch_all_by_RawContig_constraint
Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor::fetch_all_by_RawContig
Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor::fetch_all_by_RawContig_and_score
Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor::remove_by_RawContig
Bio::EnsEMBL::DBSQL::DBAdaptor::db_handle
Bio::EnsEMBL::DBSQL::DBAdaptor::port
Bio::EnsEMBL::DBSQL::DBAdaptor::driver
Bio::EnsEMBL::DBSQL::DBAdaptor::password
Bio::EnsEMBL::DBSQL::DBAdaptor::username
Bio::EnsEMBL::DBSQL::DBAdaptor::host
Bio::EnsEMBL::DBSQL::DBAdaptor::reconnect_when_lost
Bio::EnsEMBL::DBSQL::DBAdaptor::disconnect_when_inactive
Bio::EnsEMBL::DBSQL::DBAdaptor::dbname
Bio::EnsEMBL::DBSQL::DBAdaptor::prepare
Bio::EnsEMBL::DBSQL::DBAdaptor::list_supported_assemblies
Bio::EnsEMBL::DBSQL::DBAdaptor::assembly_type
Bio::EnsEMBL::DBSQL::DBAdaptor::db
Bio::EnsEMBL::DBSQL::DBAdaptor::source
Bio::EnsEMBL::DBSQL::DBConnection::group
Bio::EnsEMBL::DBSQL::DBConnection::species
Bio::EnsEMBL::DBSQL::DBEntryAdaptor::geneids_by_extids
Bio::EnsEMBL::DBSQL::DBEntryAdaptor::translationids_by_extids
Bio::EnsEMBL::DBSQL::DBEntryAdaptor::transcriptids_by_extids
Bio::EnsEMBL::DBSQL::DataFileAdaptor::DataFile_to_extension
Bio::EnsEMBL::DBSQL::ExonAdaptor::get_stable_entry_info
Bio::EnsEMBL::DBSQL::ExonAdaptor::fetch_all_by_gene_id
Bio::EnsEMBL::DBSQL::GeneAdaptor::fetch_nearest_Gene_by_Feature
Bio::EnsEMBL::DBSQL::GeneAdaptor::fetch_by_maximum_DBLink
Bio::EnsEMBL::DBSQL::GeneAdaptor::get_display_xref
Bio::EnsEMBL::DBSQL::GeneAdaptor::get_description
Bio::EnsEMBL::DBSQL::GeneAdaptor::fetch_all_by_DBEntry
Bio::EnsEMBL::DBSQL::GeneAdaptor::get_stable_entry_info
Bio::EnsEMBL::DBSQL::GeneAdaptor::fetch_by_Peptide_id
Bio::EnsEMBL::DBSQL::MetaContainer::get_Species


Removed in e84
Bio::EnsEMBL::DBSQL::CoordSystemAdaptor::_fetch_by_attrib
Bio::EnsEMBL::DBSQL::CoordSystemAdaptor::_fetch_all_by_attrib
Bio::EnsEMBL::DBSQL::MetaContainer::get_short_name
