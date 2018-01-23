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

# patch_54_55_e.sql
#
# title: Add column 'is_constitutive' to the exon table.
#
# description:
# The 'is_constitutive' column in the exon table will be set to true
# for exons that are constitutive.  This is done by a script in
# 'misc-scripts'.

ALTER TABLE exon ADD COLUMN is_constitutive BOOLEAN NOT NULL DEFAULT 0;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
VALUES (NULL, 'patch', 'patch_54_55_e.sql|add_is_constitutive_column');
