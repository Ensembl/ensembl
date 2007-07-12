CREATE TABLE `external_db` (
  `external_db_id` int(11) NOT NULL default '0',
  `db_name` varchar(27) collate latin1_bin NOT NULL default '',
  `db_release` varchar(40) collate latin1_bin NOT NULL default '',
  `status` enum('KNOWNXREF','KNOWN','XREF','PRED','ORTH','PSEUDO') collate latin1_bin NOT NULL default 'KNOWNXREF',
  `dbprimary_acc_linkable` tinyint(1) NOT NULL default '1',
  `display_label_linkable` tinyint(1) NOT NULL default '0',
  `priority` int(11) NOT NULL default '0',
  `db_display_name` varchar(255) collate latin1_bin default NULL,
  `type`  ENUM('ARRAY', 'ALT_TRANS', 'MISC', 'LIT'),
  `secondary_db_name`          VARCHAR(255),
  `secondary_db_table`         VARCHAR(255),

  PRIMARY KEY  (`external_db_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

