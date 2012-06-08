CREATE TABLE intron_supporting_evidence (
intron_supporting_evidence_id INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
analysis_id                   SMALLINT UNSIGNED NOT NULL,
seq_region_id                 INT(10) UNSIGNED NOT NULL,
seq_region_start              INT(10) UNSIGNED NOT NULL,
seq_region_end                INT(10) UNSIGNED NOT NULL,
seq_region_strand             TINYINT(2) NOT NULL,
hit_name                      VARCHAR(100) NOT NULL,
score                         DECIMAL(10,3),
score_type                    ENUM('NONE', 'DEPTH') DEFAULT 'NONE',
is_splice_canonical           BOOLEAN NOT NULL DEFAULT 0,

PRIMARY KEY (intron_supporting_evidence_id),

UNIQUE KEY (analysis_id, seq_region_id, seq_region_start, seq_region_end, seq_region_strand, hit_name)
) COLLATE=latin1_swedish_ci ENGINE=MyISAM;