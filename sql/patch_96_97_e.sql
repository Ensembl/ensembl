-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2020] EMBL-European Bioinformatics Institute
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

# patch_96_97_e.sql
#
# Title: Add rnaproduct as accepted type for stable_id_event
#
# Description:
#   This is so that Ensembl stable-id mapping can support mature RNA products, e.g. microRNA

ALTER TABLE stable_id_event MODIFY COLUMN type enum('gene', 'transcript', 'translation', 'rnaproduct') NOT NULL;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_96_97_e.sql|add_stable_id_event_type_rnaproduct');
