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

# patch_74_75_e.sql
#
# title: Add unique constraint on attrib tables
#
# description:
# For all attrib related tables, attrib_type_id and value should be unique for a given object

ALTER TABLE seq_region_attrib ADD UNIQUE KEY region_attribx (seq_region_id, attrib_type_id, value(500));
ALTER TABLE gene_attrib ADD UNIQUE KEY gene_attribx (gene_id, attrib_type_id, value(500));
ALTER TABLE transcript_attrib ADD UNIQUE KEY transcript_attribx (transcript_id, attrib_type_id, value(500));
ALTER TABLE translation_attrib ADD UNIQUE KEY translation_attribx (translation_id, attrib_type_id, value(500));
ALTER TABLE misc_attrib ADD UNIQUE KEY misc_attribx (misc_feature_id, attrib_type_id, value(500));


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_74_75_e.sql|unique_attrib_key');

 
