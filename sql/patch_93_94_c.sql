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

# patch_93_94_c.sql
#
# Title: Default align_type is NULL
#
# Description:
#   Ensure default align_type in align_feature tables is 'ensembl'

ALTER TABLE protein_feature MODIFY COLUMN align_type ENUM('ensembl', 'cigar', 'cigarplus', 'vulgar', 'mdtag') DEFAULT NULL;
UPDATE protein_feature SET align_type = NULL WHERE align_type = '';

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_93_94_c.sql|default_aln_type');
