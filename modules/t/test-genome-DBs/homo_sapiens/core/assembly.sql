CREATE TABLE `assembly` (
  `asm_seq_region_id` int(10) unsigned NOT NULL default '0',
  `cmp_seq_region_id` int(10) unsigned NOT NULL default '0',
  `asm_start` int(10) NOT NULL default '0',
  `asm_end` int(10) NOT NULL default '0',
  `cmp_start` int(10) NOT NULL default '0',
  `cmp_end` int(10) NOT NULL default '0',
  `ori` tinyint(4) NOT NULL default '0',
  KEY `cmp_seq_region_id` (`cmp_seq_region_id`),
  KEY `asm_seq_region_id` (`asm_seq_region_id`,`asm_start`),
  UNIQUE KEY all_idx (asm_seq_region_id, cmp_seq_region_id, asm_start, asm_end, cmp_start, cmp_end, ori)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

