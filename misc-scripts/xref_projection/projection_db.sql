# Table structure for projection info database

CREATE TABLE projections (

  db_release					INT NOT null,
  timestamp					DATETIME,
  from_db					VARCHAR(255),
  from_species_latin				VARCHAR(255),
  from_species_common				VARCHAR(255),
  to_db					        VARCHAR(255),
  to_species_latin				VARCHAR(255),
  to_species_common				VARCHAR(255)
  
);
