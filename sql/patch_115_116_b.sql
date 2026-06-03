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

# patch_115_116_b.sql
#
# Title: Added indices to transcript and assembly
#
# Description:
#   Ensure meta_value is not null
ALTER TABLE assembly ADD INDEX asm_overlap_end_idx (asm_seq_region_id, asm_end, asm_start);
ALTER TABLE transcript ADD INDEX seq_region_current_start_idx (seq_region_id, is_current, seq_region_start, transcript_id, gene_id);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_115_116_b.sql|Added indices to transcript and assembly');
