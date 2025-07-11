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

# patch_50_51_i.sql
#
# title: Meta value binary
#
# description:
# Add BINARY flag to meta.meta_value to make it case sensitive when doing UNIQUE comparisions (especially for Gallus/gallus)

# Drop UNIQUE index first
ALTER TABLE meta DROP INDEX species_key_value_idx;

ALTER TABLE meta CHANGE COLUMN meta_value meta_value VARCHAR(255) BINARY NOT NULL;

# Redo index
ALTER TABLE meta ADD UNIQUE INDEX species_key_value_idx (species_id, meta_key, meta_value);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_50_51_i.sql|meta_value_binary');


