CREATE TABLE ditag (

       ditag_id INT(10) NOT NULL auto_increment,
       name VARCHAR(30),
       type VARCHAR(30),
       tag_count smallint(6) default 1,
       sequence TEXT,

       PRIMARY KEY ( ditag_id )
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
