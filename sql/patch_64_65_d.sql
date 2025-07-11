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

# patch_62_63_d.sql
#
# Title: Add new Xref info_type
#
# Description:
# We can now do checksum mapping in the Xref pipeline so Xrefs need to reflect
# this

ALTER TABLE xref modify COLUMN info_type ENUM( 
                                    'PROJECTION', 'MISC', 'DEPENDENT',
                                    'DIRECT', 'SEQUENCE_MATCH',
                                    'INFERRED_PAIR', 'PROBE',
                                    'UNMAPPED', 'COORDINATE_OVERLAP',
                                    'CHECKSUM');

# Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_64_65_d.sql|add_checksum_info_type');
