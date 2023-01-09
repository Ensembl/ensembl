-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2022] EMBL-European Bioinformatics Institute
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

# patch_109_110_b.sql
#
# Title: Add 'IS_PAR' relationship
#
# Description:
#   Add 'IS_PAR' relationship to link X- and Y-PAR genes via alt allele data

ALTER TABLE alt_allele_attrib MODIFY attrib ENUM('IS_REPRESENTATIVE','IS_MOST_COMMON_ALLELE','IN_CORRECTED_ASSEMBLY','HAS_CODING_POTENTIAL','IN_ARTIFICIALLY_DUPLICATED_ASSEMBLY','IN_SYNTENIC_REGION','HAS_SAME_UNDERLYING_DNA_SEQUENCE','IN_BROKEN_ASSEMBLY_REGION','IS_VALID_ALTERNATE','SAME_AS_REPRESENTATIVE','SAME_AS_ANOTHER_ALLELE','MANUALLY_ASSIGNED','AUTOMATICALLY_ASSIGNED', 'IS_PAR');

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_109_110_b.sql|Add IS_PAR relationship to link X- and Y-PAR genes');
