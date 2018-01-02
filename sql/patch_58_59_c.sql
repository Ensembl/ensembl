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

# patch_58_59_c.sql
#
# Title:
#   Replace splicing_event.type with splicing_event.attrib_type_id
#
# Description:
#   The 'type' enumeration in the splicing_event table is too terse.
#   Replace it with a reference to a proper attrib_type.

# Add the new column.
ALTER TABLE splicing_event
  ADD COLUMN attrib_type_id SMALLINT(5) UNSIGNED NOT NULL DEFAULT 0
    AFTER seq_region_strand;

# Update the new column using info from the old column and the
# attrib_type table.
UPDATE  splicing_event se
  JOIN  attrib_type at ON (se.`type` = at.code)
SET     se.attrib_type_id = at.attrib_type_id
WHERE   at.attrib_type_id BETWEEN 300 AND 311;

# Drop the old column.
ALTER TABLE splicing_event DROP COLUMN `type`;

# Insert patch identifier.
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch',
    'patch_58_59_c.sql|splicing_event_attrib_type_id');
