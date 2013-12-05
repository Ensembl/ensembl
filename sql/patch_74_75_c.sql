# patch_74_75_c.sql
#
# title: Add genome_statistics table
#
# description:
# Addition of a new table, genome_statistics, to store genome related statistics

CREATE TABLE genome_statistics (

  genome_statistics_id     INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  statistic                VARCHAR(128) NOT NULL,
  value                    INT(10) UNSIGNED DEFAULT '0' NOT NULL,
  species_id               INT UNSIGNED DEFAULT 1,
  attrib_type_id           INT(10) UNSIGNED  DEFAULT NULL,
  timestamp                DATETIME NOT NULL DEFAULT '0000-00-00 00:00:00',

  PRIMARY KEY (genome_statistics_id),
  UNIQUE KEY stats_uniq(statistic, attrib_type_id, species_id),
  KEY stats_idx (statistic, attrib_type_id, species_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_74_75_c.sql|add_genome_statistics');

 
