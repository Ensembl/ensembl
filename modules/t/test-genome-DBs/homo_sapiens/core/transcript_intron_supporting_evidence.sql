CREATE TABLE transcript_intron_supporting_evidence (
transcript_id                 INT(10) UNSIGNED NOT NULL,
intron_supporting_evidence_id INT(10) UNSIGNED NOT NULL,
previous_exon_id              INT(10) UNSIGNED NOT NULL,
next_exon_id                  INT(10) UNSIGNED NOT NULL,
PRIMARY KEY (intron_supporting_evidence_id, transcript_id)
) COLLATE=latin1_swedish_ci ENGINE=MyISAM;