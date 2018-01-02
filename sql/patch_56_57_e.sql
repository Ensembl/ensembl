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

# patch_56_57_e.sql
#
# Title: Add canonical translations
#
# Description:
# Add a reference to the translation table in the transcript table to
# let the API to know which translation (out of possibly several) is the
# canonical one.
#
# This patch adds a 'canonical_translation_id' field (with an index)
# to the 'transcript' table and initiates it with the corresponding
# 'translation.transcript_id'.

ALTER TABLE transcript
ADD COLUMN canonical_translation_id INT(10) UNSIGNED,
ADD UNIQUE INDEX canonical_translation_idx (canonical_translation_id);

-- Initiate the new transcript.canonical_translation_id field with
-- translation.transcript_id.
UPDATE  transcript ts, translation tl
SET     ts.canonical_translation_id = tl.translation_id
WHERE   ts.transcript_id = tl.transcript_id;

OPTIMIZE TABLE transcript;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
VALUES (NULL, 'patch', 'patch_56_57_e.sql|canonical_translations');
