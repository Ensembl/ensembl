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

# patch_54_55_g.sql
#
# title: Change analysis_description.display_label to "NOT NULL".
#
# description:
# Display labels are generally supposed to be non-NULL and non-empty.
# The display_label field in the analysis_description table allows for
# NULLs.  This patch fixes this.

ALTER TABLE analysis_description
  MODIFY COLUMN display_label VARCHAR(255) NOT NULL;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
VALUES (NULL, 'patch',
  'patch_54_55_g.sql|analysis_description.display_label_NOT_NULL');
