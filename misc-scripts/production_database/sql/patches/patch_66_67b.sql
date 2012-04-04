-- Patching in support for numerous different names for trinomials

-- Supporting aliases
CREATE TABLE species_alias (
  species_alias_id  INTEGER UNSIGNED NOT NULL AUTO_INCREMENT,
  species_id        INTEGER UNSIGNED NOT NULL,
  alias             varchar(255) NOT NULL,
  is_current        BOOLEAN NOT NULL DEFAULT true,
  created_by    INTEGER,
  created_at    DATETIME,
  modified_by   INTEGER,
  modified_at   DATETIME,
  PRIMARY KEY (species_alias_id),
  UNIQUE INDEX (alias, is_current),
  INDEX sa_speciesid_idx (species_id)
);


-- Supporting production name
ALTER TABLE species
ADD COLUMN production_name VARCHAR(255) NOT NULL
AFTER web_name;

-- Supporting scientific name
ALTER TABLE species
ADD COLUMN scientific_name VARCHAR(255) NOT NULL
AFTER web_name;

-- Supporting URL name
ALTER TABLE species
ADD COLUMN url_name VARCHAR(255) NOT NULL
AFTER production_name;