CREATE TABLE mapping_set (

        mapping_set_id  INT(10) UNSIGNED NOT NULL,
        schema_build    VARCHAR(20) NOT NULL,

        PRIMARY KEY(schema_build)

) ENGINE=MyISAM DEFAULT CHARSET=latin1;
