-- See the NOTICE file distributed with this work for additional information
-- regarding copyright ownership.
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

# patch_55_56_d.sql
#
# Title: Add an index to the splicing_event_feature table
#
# Description:
# With an index on transcript_id in splicing_event_feature, the
# generation of biomarts will be sped up.

ALTER TABLE splicing_event_feature
  ADD INDEX transcript_idx (transcript_id);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
VALUES (NULL, 'patch', 'patch_55_56_d.sql|add_index_to_splicing_event_feature');
