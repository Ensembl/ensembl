DROP INDEX stable_id_idx on exon_stable_id;
CREATE INDEX  stable_id_idx on exon_stable_id(stable_id, version);
DROP INDEX stable_id_idx on gene_stable_id;
CREATE INDEX  stable_id_idx on gene_stable_id(stable_id, version);
DROP INDEX stable_id_idx on transcript_stable_id;
CREATE INDEX  stable_id_idx on transcript_stable_id(stable_id, version);
DROP INDEX stable_id_idx on translation_stable_id;
CREATE INDEX  stable_id_idx on translation_stable_id(stable_id, version);
