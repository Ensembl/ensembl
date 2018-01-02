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

# patch_51_52_b.sql
#
# title: Widen some text columns
#
# description:
# Change analysis.parameters to text, xref.description to text, coord_system.version to varchar(255), DNA.sequence to LONGTEXT

ALTER TABLE analysis MODIFY parameters TEXT;

ALTER TABLE xref MODIFY description TEXT;

ALTER TABLE coord_system MODIFY version VARCHAR(255) DEFAULT NULL;

ALTER TABLE dna MODIFY sequence LONGTEXT NOT NULL;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_51_52_b.sql|widen_columns');


