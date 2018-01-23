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

# patch_60_61_c.sql
#
# Title:
#   Modify the indexes of the object_xref tables.
#
# Description:
# The object_xref_id field of the object_xref table was found to for
# some reason not be a primary key.  This field is always unique, so why
# not add it as the primary key?  Also re-jig the other indexes while we
# are at it.

ALTER TABLE object_xref
  ADD PRIMARY KEY (object_xref_id);

# Also notice how the existing oxref_idx index now duplicates the
# object_xref_id index (although it is not "UNIQUE", it really is).
# Removing object_xref_id out of this index, it now is a wider version
# of the xref_idx index.  After some testing, we realize we can drop the
# oxref_idx index completely (the API never uses it), reorder and rename
# the unique object_type_idx index and replace the old xref_idx with an
# index covering iensembl_object_type and ensembl_id.  The paricular
# ordering that we have chosen seems to be the most optimal for the use
# cases that we have in the API.

ALTER TABLE object_xref
  DROP INDEX object_type_idx,
  DROP INDEX oxref_idx,
  DROP INDEX xref_idx,
  ADD UNIQUE INDEX
    xref_idx (xref_id, ensembl_object_type, ensembl_id, analysis_id),
  ADD INDEX ensembl_idx (ensembl_object_type, ensembl_id);

OPTIMIZE TABLE object_xref;

# Insert patch identifier.
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_60_61_c.sql|rejig_object_xref_indexes');

