# patch_72_73_c.sql
#
# title: Add Marker as accepted ensembl_object_type for object_xref
#
# description:
# Add missing GO_XREF enums that are missing

ALTER TABLE object_xref MODIFY COLUMN ensembl_object_type enum('RawContig', 'Transcript', 'Gene', 'Translation', 'Operon', 'OperonTranscript', 'Marker') NOT NULL;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_72_73_c.sql|add_object_type_marker');

 
