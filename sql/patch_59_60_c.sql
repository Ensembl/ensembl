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

# patch_59_60_c.sql
#
# Title:
#   A patch to fix a couple of inconsistencies in the schema.
#
# Description:
#   QC turned up issues with the signedness of a number of fields in
#   the Ensembl Core schema.  This patch fixes these.  The fields
#   are: karyotype.seq_region_start, karyotype.seq_region_end and
#   seq_region.length (should all be UNSIGNED).

# Make the the 'seq_region_start' and 'seq_region_end' fields of the
# 'karyotype' table UNSIGNED (like they are everywhere else).
ALTER TABLE karyotype
  MODIFY COLUMN seq_region_start INT(10) UNSIGNED NOT NULL,
  MODIFY COLUMN seq_region_end   INT(10) UNSIGNED NOT NULL;

# Make 'seq_region.length' UNSIGNED (we do not like negative lengths).
ALTER TABLE seq_region
  MODIFY COLUMN length INT(10) UNSIGNED NOT NULL;

# Insert patch identifier.
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_59_60_c.sql|QC_fixes');
