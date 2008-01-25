# patch_48_49_e.sql
#
# ensembl_object_type_not_null
#
# description:
# ensembl_object_type not null

ALTER TABLE object_xref CHANGE COLUMN ensembl_object_type
  ensembl_object_type ENUM('RawContig', 'Transcript', 'Gene',
                                   'Translation') NOT NULL;

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_48_49_e.sql|ensembl_object_type_not_null');



