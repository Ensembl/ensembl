CREATE TABLE `unconventional_transcript_association` (

       `transcript_id`    INT(10) UNSIGNED NOT NULL,
       `gene_id`          INT(10) UNSIGNED NOT NULL,
       `interaction_type` ENUM("antisense","sense_intronic","sense_overlaping_exonic","chimeric_sense_exonic"),

       KEY (transcript_id),
       KEY (gene_id)

) ENGINE=MyISAM DEFAULT CHARSET=latin1;
