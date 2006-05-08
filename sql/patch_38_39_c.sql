# patch_38_39_c
#
# title: multi-version objects
#
# description:
# this patch will allow the storage of multiple versions of genes,
# transcripts and exons in a single database

# add column 'is_current' boolean default 1

ALTER TABLE gene ADD COLUMN is_current BOOLEAN DEFAULT 1;
ALTER TABLE transcript ADD COLUMN is_current BOOLEAN DEFAULT 1;
ALTER TABLE exon ADD COLUMN is_current BOOLEAN DEFAULT 1;

# change UNIQUE KEY 'stable_id' to normal KEY in stable_id tables

ALTER TABLE gene_stable_id DROP INDEX stable_id;
ALTER TABLE gene_stable_id ADD INDEX stable_id_idx (stable_id, version);

ALTER TABLE transcript_stable_id DROP INDEX stable_id;
ALTER TABLE transcript_stable_id ADD INDEX stable_id_idx (stable_id, version);

ALTER TABLE exon_stable_id DROP INDEX stable_id;
ALTER TABLE exon_stable_id ADD INDEX stable_id_idx (stable_id, version);

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_38_39_c.sql|multiversion_objects');

