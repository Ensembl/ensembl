# Patch to convert release 37 Ensembl schema to release 38

UPDATE meta set meta_value="38" where meta_key="schema_version";

# Add info_type & info_text columns to xref

ALTER TABLE xref ADD COLUMN info_type ENUM('PROJECTION', 'MISC');
ALTER TABLE xref ADD COLUMN info_text VARCHAR(255);
ALTER TABLE xref ADD INDEX info_type_idx (info_type);

# Change name of release column in external_db since release is a reserved word
ALTER TABLE external_db CHANGE COLUMN release db_release VARCHAR(40) NOT NULL;

# Add the two new Unmapped Object tables:-


################################################################################
#
# Table structure for table 'unmapped_object'
#
# Describes why a particular external entity was not mapped to an ensembl one.

CREATE TABLE unmapped_object (

  unmapped_object_id    INT UNSIGNED NOT NULL AUTO_INCREMENT,
  type                  ENUM('xref', 'cDNA', 'Marker') NOT NULL,
  analysis_id           INT(10) UNSIGNED NOT NULL,
  external_db_id        INT NOT NULL,
  identifier            VARCHAR(255) NOT NULL,
  unmapped_reason_id    SMALLINT(5) UNSIGNED NOT NULL,
  query_score           DOUBLE,
  target_score          DOUBLE,
  ensembl_id            INT(10) unsigned default '0',
  ensembl_object_type   ENUM('RawContig','Transcript','Gene','Translation') collate latin1_bin default 'RawContig',
  PRIMARY KEY            ( unmapped_object_id ),
  KEY                    id_idx( identifier ),
  KEY                    anal_idx( analysis_id ),
  KEY                    anal_exdb_idx( analysis_id, external_db_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;

################################################################################
#
# Table structure for table 'unmapped_reason'
#
# Describes the reason why a mapping failed.

CREATE TABLE unmapped_reason (

  unmapped_reason_id     SMALLINT(5) UNSIGNED NOT NULL AUTO_INCREMENT,
  summary_description    VARCHAR(255),
  full_description       VARCHAR(255),

  PRIMARY KEY ( unmapped_reason_id )

) COLLATE=latin1_swedish_ci TYPE=MyISAM;
