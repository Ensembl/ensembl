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

# patch_38_39_a
#
# title: status enum
#
# description:
# this patch adds a new status to the enumeration in gene & transcript

# Add a new status to the enumeration in gene & transcript

ALTER TABLE gene CHANGE COLUMN status status enum('KNOWN','NOVEL','PUTATIVE','PREDICTED','KNOWN_BY_PROJECTION');

ALTER TABLE transcript CHANGE COLUMN status status enum('KNOWN','NOVEL','PUTATIVE','PREDICTED','KNOWN_BY_PROJECTION');

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_38_39_a.sql|status_enum');

