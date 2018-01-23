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

# patch_52_53_c.sql
#
# title: identity xref rename
#
# description:
# Rename some columns in identity xref to make them clearer
#
# query_identity    ->  xref_identity
# target_identity   ->  ensembl_identity
# hit_start         ->  xref_start
# hit_end           ->  xref_end
# translation_start ->  ensembl_start
# translation_end   ->  ensembl_end

ALTER TABLE identity_xref CHANGE COLUMN query_identity xref_identity INT(5);
ALTER TABLE identity_xref CHANGE COLUMN target_identity ensembl_identity INT(5);
ALTER TABLE identity_xref CHANGE COLUMN hit_start xref_start INT;
ALTER TABLE identity_xref CHANGE COLUMN hit_end xref_end INT;
ALTER TABLE identity_xref CHANGE COLUMN translation_start ensembl_start INT;
ALTER TABLE identity_xref CHANGE COLUMN translation_end ensembl_end INT;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_52_53_c.sql|identity_xref_rename');


