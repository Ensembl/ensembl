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

# patch_53_54_c.sql
#
# title: Move analysis_id from identity_xref to object_xref
#
# description:
# Add analysis_xref column to object_xref, copy values from identity_xref, remove column from identity_xref.
# Will allow all object_xref relationships to have an analysis, not just those from sequence matching.

ALTER TABLE object_xref ADD COLUMN analysis_id SMALLINT UNSIGNED DEFAULT NULL;

UPDATE object_xref ox, identity_xref ix SET ox.analysis_id=ix.analysis_id WHERE ox.object_xref_id=ix.object_xref_id;

ALTER TABLE object_xref ADD KEY analysis_idx (analysis_id);

ALTER TABLE identity_xref DROP KEY analysis_idx;

ALTER TABLE identity_xref DROP COLUMN analysis_id;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_53_54_c.sql|identity_object_analysis_move');


