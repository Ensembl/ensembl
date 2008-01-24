# patch_48_49_b.sql
#
# title: new canonical transcript column n gene table
#
# description:
# Add canonical_transcript column to gene table

ALTER TABLE gene ADD COLUMN canonical_transcript INT(10) UNSIGNED;

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_48_49_b.sql|new_canonical_transcript_column');
