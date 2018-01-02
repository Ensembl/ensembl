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

# patch_79_80_b.sql
#
# Title: Extend column width for xref.dbprimary_acc
#
# Description:
#   Column width extended from VARCHAR(40) to VARCHAR(50) for xref.dbprimary_acc

ALTER TABLE xref MODIFY COLUMN dbprimary_acc VARCHAR(50) NOT NULL;

# Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_79_80_b.sql|xref_dbprimary_acc_longer');
