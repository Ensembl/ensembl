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

# patch_63_64_d.sql
#
# Title:
#   change to linkage_type column in ontology_xref 
#
# Description:
#    Change the column definition to VARCHAR(3) instead of providing
#    a list of allowable values. Values will be tested in healthchecks.

ALTER TABLE ontology_xref DROP INDEX object_source_type_idx;

ALTER TABLE ontology_xref
  MODIFY linkage_type VARCHAR(3) DEFAULT NULL;

ALTER TABLE ontology_xref 
  ADD UNIQUE INDEX object_source_type_idx (object_xref_id, source_xref_id, linkage_type);

# Insert patch identifier.
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch',
    'patch_63_64_d.sql|linkage_type change in ontology_xref');
