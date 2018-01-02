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

# patch_39_40_e
#
# title: is_current not null
#
# description:
# this patch changes the is_current column in gene, transcript and exon to be
# 'NOT NULL', to prevent loading data where this property is not set

# change column 'is_current'

ALTER TABLE gene CHANGE is_current is_current BOOLEAN NOT NULL DEFAULT 1;
ALTER TABLE transcript CHANGE is_current is_current BOOLEAN NOT NULL DEFAULT 1;
ALTER TABLE exon CHANGE is_current is_current BOOLEAN NOT NULL DEFAULT 1;

# set is_current to 1
UPDATE gene set is_current = 1;
UPDATE transcript set is_current = 1;
UPDATE exon set is_current = 1;

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_39_40_e.sql|is_current_not_null');

