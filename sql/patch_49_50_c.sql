# patch_49_50_c.sql
#
# title: Modify canonical transcript support
#
# description:
# Change gene.canonical_transcript to canonical_transcript_id, and add canonical_annotation, DEFAULT NULL

ALTER TABLE gene CHANGE COLUMN canonical_transcript canonical_transcript_id INT(10) UNSIGNED NOT NULL;
ALTER TABLE gene ADD COLUMN canonical_annotation VARCHAR(255) DEFAULT NULL;

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_49_50_c.sql|canonical_transcript');


