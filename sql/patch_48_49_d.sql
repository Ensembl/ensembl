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

# patch_48_49_d.sql
#
# Title: Add new info_type to xref table
#
# Description:
# Add the ENUM 'COORDINATE_OVERLAP' to the xref.info_type column.

ALTER TABLE xref CHANGE COLUMN info_type
  info_type ENUM(
    'PROJECTION', 'MISC', 'DEPENDENT', 'DIRECT', 'SEQUENCE_MATCH',
    'INFERRED_PAIR', 'PROBE', 'UNMAPPED', 'COORDINATE_OVERLAP'
  );

# Patch identifier
INSERT INTO meta (meta_key, meta_value)
VALUES ('patch', 'patch_48_49_d.sql|new_info_type_enum');
