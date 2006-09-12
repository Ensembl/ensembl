# patch_40_41_f
#
# title: attrib tavle indices
#
# description:
# Add extra indices to attrib tables to improve performance.

ALTER TABLE misc_attrib ADD INDEX val_only_idx (value(40));
ALTER TABLE seq_region_attrib ADD INDEX val_only_idx (value(40));
ALTER TABLE gene_attrib ADD INDEX val_only_idx (value(40));
ALTER TABLE transcript_attrib ADD INDEX val_only_idx (value(40));
ALTER TABLE translation_attrib ADD INDEX val_only_idx (value(40));

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_40_41_f.sql|attrib_indices');


