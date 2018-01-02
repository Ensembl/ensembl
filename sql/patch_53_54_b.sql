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

# patch_53_54_b.sql
#
# title: Widen some text columns
#
# description:
# Change oligo_probe.name to 40 characters, external_db.db_name to 100 chars, analysis.logic_name to 128 characters

ALTER TABLE oligo_probe MODIFY name VARCHAR(40);

ALTER TABLE external_db MODIFY db_name VARCHAR(100) NOT NULL;

ALTER TABLE analysis MODIFY logic_name VARCHAR(128) NOT NULL;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_53_54_b.sql|widen_columns');


