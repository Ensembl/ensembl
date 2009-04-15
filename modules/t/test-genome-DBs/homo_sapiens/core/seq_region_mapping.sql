CREATE TABLE seq_region_mapping (

        external_seq_region_id  INT(10) UNSIGNED NOT NULL,
        internal_seq_region_id  INT(10) UNSIGNED NOT NULL,
        mapping_set_id          INT(10) UNSIGNED NOT NULL,

        KEY (mapping_set_id)

) ENGINE=MyISAM DEFAULT CHARSET=latin1;

