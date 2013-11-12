# patch_74_75_b.sql
#
# Title: Add source to transcript
#
# Description: Column allowing an optional source to be added to a transcript

ALTER TABLE transcript ADD COLUMN source VARCHAR(20) NOT NULL;

# Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_74_75_b.sql|transcript_source');


