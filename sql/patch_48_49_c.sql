# patch_48_49_c.sql
#
# title: regulatory_support_removal
#
# description:
# regulatory tables to be removed from database (now done by func gen)

DELETE object_xref FROM object_xref where ensembl_object_type = "regulatory_factor";
DELETE object_xref FROM object_xref where ensembl_object_type = "regulatory_feature";

ALTER TABLE object_xref CHANGE COLUMN ensembl_object_type
  ensembl_object_type ENUM('RawContig', 'Transcript', 'Gene',
                                   'Translation');

DROP TABLE regulatory_factor;
DROP TABLE regulatory_factor_coding;
DROP TABLE regulatory_feature;
DROP TABLE regulatory_feature_object;
DROP TABLE regulatory_search_region;

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_48_49_c.sql|regulatory_support_removal');



