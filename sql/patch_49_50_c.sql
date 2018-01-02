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

# patch_49_50_c.sql
#
# title: Modify canonical transcript support
#
# description:
# Change gene.canonical_transcript to canonical_transcript_id, and add canonical_annotation, DEFAULT NULL

ALTER TABLE gene CHANGE COLUMN canonical_transcript canonical_transcript_id INT(10) UNSIGNED NOT NULL;
ALTER TABLE gene ADD COLUMN canonical_annotation VARCHAR(255) DEFAULT NULL;

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_49_50_c.sql|canonical_transcript');


