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

# patch_38_39_c
#
# title: multi-version objects
#
# description:
# this patch will allow the storage of multiple versions of genes,
# transcripts and exons in a single database

# add column 'is_current' boolean default 1

ALTER TABLE gene ADD COLUMN is_current BOOLEAN DEFAULT 1;
ALTER TABLE transcript ADD COLUMN is_current BOOLEAN DEFAULT 1;
ALTER TABLE exon ADD COLUMN is_current BOOLEAN DEFAULT 1;

# change UNIQUE KEY 'stable_id' to normal KEY in stable_id tables

ALTER TABLE gene_stable_id DROP INDEX stable_id;
ALTER TABLE gene_stable_id ADD INDEX stable_id_idx (stable_id, version);

ALTER TABLE transcript_stable_id DROP INDEX stable_id;
ALTER TABLE transcript_stable_id ADD INDEX stable_id_idx (stable_id, version);

ALTER TABLE exon_stable_id DROP INDEX stable_id;
ALTER TABLE exon_stable_id ADD INDEX stable_id_idx (stable_id, version);

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_38_39_c.sql|multiversion_objects');

