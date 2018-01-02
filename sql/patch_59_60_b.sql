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

# patch_59_60_b.sql
#
# Title:
#   Rename 'go_xref' table to 'ontology_xref'.
#
# Description:
#   Rename the 'go_xref' table to make its use more generic.

# Rename the table, and swap the source_xref_id and linkage_type fields.
ALTER TABLE go_xref
  RENAME TO ontology_xref,
  MODIFY COLUMN source_xref_id INT(10) UNSIGNED DEFAULT NULL
    AFTER object_xref_id;

# Optimize the table, because indexes may be out of whack
OPTIMIZE TABLE ontology_xref;

# Insert patch identifier.
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_59_60_b.sql|rename_go_xref_table');
