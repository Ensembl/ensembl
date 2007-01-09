# patch_42_43_e
#
# title: unmapped_object.parent
#
# description:
# Add an index on gene_archive.peptide_archive_id

ALTER TABLE gene_archive ADD INDEX peptide_archive_id_idx (peptide_archive_id);

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_42_43_e.sql|gene_archive.peptide_archive_id.index');


