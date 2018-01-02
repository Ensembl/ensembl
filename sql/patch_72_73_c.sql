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

# patch_72_73_c.sql
#
# title: Add Marker as accepted ensembl_object_type for object_xref
#
# description:
# Add missing GO_XREF enums that are missing

ALTER TABLE object_xref MODIFY COLUMN ensembl_object_type enum('RawContig', 'Transcript', 'Gene', 'Translation', 'Operon', 'OperonTranscript', 'Marker') NOT NULL;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_72_73_c.sql|add_object_type_marker');

 
