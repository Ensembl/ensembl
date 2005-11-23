# Patch to convert release 35 Ensembl schema to release 36

UPDATE meta set meta_value="36" where meta_key="schema_version";

# New table (so no data conversion required) for storing regulatory search regions

CREATE TABLE regulatory_search_region (

  regulatory_search_region_id  INT NOT NULL auto_increment,
  name                         VARCHAR(255) NOT NULL,
  seq_region_id                INT NOT NULL,                 # FK refs seq_region
  seq_region_start             INT NOT NULL,
  seq_region_end               INT NOT NULL,
  seq_region_strand            TINYINT NOT NULL,
  ensembl_object_type          ENUM( 'Transcript', 'Translation', 'Gene') NOT NULL,
  ensembl_object_id            INT,           # FK to gene/transcript/translation
  analysis_id                  INT NOT NULL,  # FK to analysis

  PRIMARY KEY (regulatory_search_region_id),
  KEY rsr_idx (regulatory_search_region_id),
  KEY ensembl_object_idx (ensembl_object_type, ensembl_object_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;

ALTER TABLE density_type ADD COLUMN region_features INT NOT NULL AFTER block_size;
ALTER TABLE density_type DROP INDEX analysis_id;
ALTER TABLE density_type ADD UNIQUE analysis_id( analysis_id, block_size, region_features);

ALTER TABLE regulatory_feature ADD INDEX seq_region_idx(seq_region_id, analysis_id, seq_region_start);
ALTER TABLE regulatory_feature ADD INDEX seq_region_idx_2(seq_region_id, seq_region_start);

ALTER TABLE regulatory_search_region ADD INDEX seq_region_idx(seq_region_id, seq_region_start);
ALTER TABLE regulatory_search_region ADD INDEX seq_region_idx_2(seq_region_id, seq_region_start);

ALTER TABLE go_xref modify linkage_type  enum('IC','IDA','IEA','IEP','IGI','IMP','IPI','ISS','NAS','ND','TAS','NR','RCA') NULL;