-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2019] EMBL-European Bioinformatics Institute
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

# patch_93_94_b.sql
#
# Title: Allow NULL analysis_id in object_xref table
#
# Description:
#   Allow analysis_id to not be set in the object_xref table, as not all entries will have an analysis

ALTER TABLE object_xref MODIFY COLUMN analysis_id SMALLINT UNSIGNED;
UPDATE object_xref SET analysis_id = NULL WHERE analysis_id = 0;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_93_94_b.sql|nullable_ox_analysis');
