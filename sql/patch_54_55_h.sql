# patch_54_55_h.sql
#
# title: Change translation-related fields in gene_archive to allow for NULLs.
#
# description:
# The gene_archive table holds deprecated genes.  The fields
# translation_stable_id, translation_version, peptide_archive_id ought
# to be NULL for genes that lack a translation.  This patch fixes these
# fields to allow for NULLs.  The patch also changes any entry with an
# empty string as translation_stable_id into NULLs for all these fields.

ALTER TABLE gene_archive
  MODIFY COLUMN translation_stable_id   VARCHAR(128),
  MODIFY COLUMN translation_version     SMALLINT,
  MODIFY COLUMN peptide_archive_id      INT(10) UNSIGNED;

-- Also fixup existing data
UPDATE gene_archive
  SET translation_stable_id = NULL,
      translation_version = NULL,
      peptide_archive_id = NULL
  WHERE translation_stable_id = '';

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
VALUES (NULL, 'patch', 'patch_54_55_h.sql|gene_archive.allow_for_NULLs');
