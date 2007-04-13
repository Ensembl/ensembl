CREATE TABLE `unmapped_reason` (

  `unmapped_reason_id`    SMALLINT(5) UNSIGNED NOT NULL AUTO_INCREMENT,
  `summary_description`   VARCHAR(255),
  `full_description`      VARCHAR(255),

  PRIMARY KEY (`unmapped_reason_id`)

) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;
