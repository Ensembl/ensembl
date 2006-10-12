# patch_41_42_b
#
# title: unconventional transcripts
#
# description:
# Add new table to support unconventional transcripts
# Also allow transcripts to not have genes (remove NOT NULL constraint on transcript.gene_id)

ALTER TABLE transcript CHANGE COLUMN gene_id gene_id INT(10) UNSIGNED;

CREATE TABLE unconventional_transcript_association (

  transcript_id    INT(10) UNSIGNED NOT NULL,
  gene_id          INT(10) UNSIGNED NOT NULL,
  interaction_type ENUM("antisense","sense_intronic","sense_overlaping_exonic","chimeric_sense_exonic"),

  KEY (transcript_id),
  KEY (gene_id)

) ENGINE=MyISAM DEFAULT CHARSET=latin1;


# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_41_42_b.sql|unconventional_transcripts');