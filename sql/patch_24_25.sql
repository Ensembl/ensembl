# Add and populate created_date and modified_date to stable_id tables
# Schema 24-25

set @today = concat( curdate(), " 12:00:00" );

ALTER TABLE gene_stable_id ADD COLUMN created_date DATETIME NOT NULL;
ALTER TABLE gene_stable_id ADD COLUMN modified_date DATETIME NOT NULL;
UPDATE gene_stable_id SET created_date=@today;
UPDATE gene_stable_id SET modified_date=@today;

ALTER TABLE exon_stable_id ADD COLUMN created_date DATETIME NOT NULL;
ALTER TABLE exon_stable_id ADD COLUMN modified_date DATETIME NOT NULL;
UPDATE exon_stable_id SET created_date=@today;
UPDATE exon_stable_id SET modified_date=@today;

ALTER TABLE transcript_stable_id ADD COLUMN created_date DATETIME NOT NULL;
ALTER TABLE transcript_stable_id ADD COLUMN modified_date DATETIME NOT NULL;
UPDATE transcript_stable_id SET created_date=@today;
UPDATE transcript_stable_id SET modified_date=@today;

ALTER TABLE translation_stable_id ADD COLUMN created_date DATETIME NOT NULL;
ALTER TABLE translation_stable_id ADD COLUMN modified_date DATETIME NOT NULL;
UPDATE translation_stable_id SET created_date=@today;
UPDATE translation_stable_id SET modified_date=@today;
