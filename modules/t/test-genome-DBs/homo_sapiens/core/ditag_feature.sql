CREATE TABLE ditag_feature (

       ditag_feature_id int(10) unsigned NOT NULL auto_increment,
       ditag_id int(10) unsigned NOT NULL default '0',
       ditag_pair_id int(10) unsigned NOT NULL default '0',
       seq_region_id int(10) unsigned NOT NULL default '0',
       seq_region_start int(10) unsigned NOT NULL default '0',
       seq_region_end int(10) unsigned NOT NULL default '0',
       seq_region_strand tinyint(1) NOT NULL default '0',
       analysis_id int(10) unsigned NOT NULL default '0',
       hit_start int(10) unsigned NOT NULL default '0',
       hit_end int(10) unsigned NOT NULL default '0',
       hit_strand tinyint(1) NOT NULL default '0',
       cigar_line text default '',
       ditag_side char default '',

       PRIMARY KEY  (ditag_feature_id),
       KEY (ditag_id),
       KEY (ditag_pair_id)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
