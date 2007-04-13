CREATE TABLE mapping_session (

  mapping_session_id 	      INT(11) NOT NULL AUTO_INCREMENT,
  old_db_name                 VARCHAR(80) NOT NULL DEFAULT '',
  new_db_name                 VARCHAR(80) NOT NULL DEFAULT '',
  old_release                 VARCHAR(5) NOT NULL DEFAULT '',
  new_release                 VARCHAR(5) NOT NULL DEFAULT '',
  old_assembly                VARCHAR(20) NOT NULL DEFAULT '',
  new_assembly                VARCHAR(20) NOT NULL DEFAULT '',
  created                     DATETIME NOT NULL,

  PRIMARY KEY (mapping_session_id)

) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

