# Add and populate created_date and modified_date to stable_id tables
# Schema 24-25

ALTER TABLE gene_stable_id ADD COLUMN created_date DATETIME NOT NULL;
ALTER TABLE gene_stable_id ADD COLUMN modified_date DATETIME NOT NULL;
UPDATE gene_stable_id SET created_date='2004-09-20 00:00:00';
UPDATE gene_stable_id SET modified_date='2004-09-20 00:00:00';

ALTER TABLE exon_stable_id ADD COLUMN created_date DATETIME NOT NULL;
ALTER TABLE exon_stable_id ADD COLUMN modified_date DATETIME NOT NULL;
UPDATE exon_stable_id SET created_date='2004-09-20 00:00:00';
UPDATE exon_stable_id SET modified_date='2004-09-20 00:00:00';

ALTER TABLE transcript_stable_id ADD COLUMN created_date DATETIME NOT NULL;
ALTER TABLE transcript_stable_id ADD COLUMN modified_date DATETIME NOT NULL;
UPDATE transcript_stable_id SET created_date='2004-09-20 00:00:00';
UPDATE transcript_stable_id SET modified_date='2004-09-20 00:00:00';

ALTER TABLE translation_stable_id ADD COLUMN created_date DATETIME NOT NULL;
ALTER TABLE translation_stable_id ADD COLUMN modified_date DATETIME NOT NULL;
UPDATE translation_stable_id SET created_date='2004-09-20 00:00:00';
UPDATE translation_stable_id SET modified_date='2004-09-20 00:00:00';
