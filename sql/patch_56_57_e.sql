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
