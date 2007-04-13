CREATE TABLE `stable_id_event` (
  `old_stable_id` varchar(128) collate latin1_bin default NULL,
  `old_version` smallint(6) default NULL,
  `new_stable_id` varchar(128) collate latin1_bin default NULL,
  `new_version` smallint(6) default NULL,
  `mapping_session_id` int(10) NOT NULL default '0',
  `type` enum('gene','transcript','translation') collate latin1_bin NOT NULL default 'gene',
  `score` FLOAT NOT NULL DEFAULT 0,
  UNIQUE KEY `uni_idx` (`mapping_session_id`,`old_stable_id`,`old_version`,`new_stable_id`,`new_version`,`type`),
  KEY `new_idx` (`new_stable_id`),
  KEY `old_idx` (`old_stable_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

