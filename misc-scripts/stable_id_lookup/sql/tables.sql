CREATE TABLE species (
  species_id		INTEGER UNSIGNED NOT NULL AUTO_INCREMENT,
  name         	     	VARCHAR(255) NOT NULL,
  taxonomy_id        	INTEGER UNSIGNED,

  PRIMARY KEY (species_id),
  UNIQUE INDEX name_idx (name)
) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


CREATE TABLE stable_id_lookup (
 stable_id   	  VARCHAR(128) NOT NULL,	      
 species_id	  INTEGER UNSIGNED NOT NULL,
 db_type          VARCHAR(255) NOT NULL,
 object_type   	  VARCHAR(255) NOT NULL,

 UNIQUE INDEX stable_id_lookup_idx (stable_id,species_id,db_type,object_type),
 KEY stable_id_db_type (stable_id,db_type,object_type),
 KEY stable_id_object_type (stable_id,object_type)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


CREATE TABLE meta (

  meta_id                     INT NOT NULL AUTO_INCREMENT,
  species_id                  INT UNSIGNED DEFAULT 1,
  meta_key                    VARCHAR(40) NOT NULL,
  meta_value                  VARCHAR(255) BINARY,

  PRIMARY   KEY (meta_id),
  UNIQUE    KEY species_key_value_idx (species_id, meta_key, meta_value),
            KEY species_value_idx (species_id, meta_value)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

