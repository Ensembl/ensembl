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



