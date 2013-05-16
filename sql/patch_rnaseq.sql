CREATE TABLE gene_expression (
  gene_expression_id          INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  gene_id                     INT(10) UNSIGNED NOT NULL,
  tissue_id                   INT(10) UNSIGNED NOT NULL,
  value                       TEXT NOT NULL,
  analysis_id                 SMALLINT UNSIGNED NOT NULL,
  value_type                  ENUM('count', 'RPKM') NOT NULL,

  PRIMARY KEY (gene_expression_id),
  UNIQUE KEY gene_expression_idx(gene_id, tissue_id, analysis_id, value_type)



) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

CREATE TABLE transcript_expression (
  transcript_expression_id          INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  transcript_id                     INT(10) UNSIGNED NOT NULL,
  tissue_id                   INT(10) UNSIGNED NOT NULL,
  value                       TEXT NOT NULL,
  analysis_id                 SMALLINT UNSIGNED NOT NULL,
  value_type                  ENUM('count', 'RPKM') NOT NULL,

  PRIMARY KEY (transcript_expression_id),
  UNIQUE KEY transcript_expression_idx(transcript_id, tissue_id, analysis_id, value_type)



) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

CREATE TABLE exon_expression (
  exon_expression_id          INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  exon_id                     INT(10) UNSIGNED NOT NULL,
  tissue_id                   INT(10) UNSIGNED NOT NULL,
  value                       TEXT NOT NULL,
  analysis_id                 SMALLINT UNSIGNED NOT NULL,
  value_type                  ENUM('count', 'RPKM') NOT NULL,

  PRIMARY KEY (exon_expression_id),
  UNIQUE KEY exon_expression_idx(exon_id, tissue_id, analysis_id, value_type)



) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

CREATE TABLE tissue (
  tissue_id                   INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  ontology                    VARCHAR(64) NOT NULL,
  name                        VARCHAR(255),
  description                 TEXT,

  PRIMARY KEY (tissue_id),
  UNIQUE KEY name_idx (name)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

