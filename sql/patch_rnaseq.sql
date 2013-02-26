CREATE TABLE gene_expression (
  gene_expression_id          INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  gene_id                     INT(10) UNSIGNED NOT NULL,
  tissue_id                   INT(10) UNSIGNED NOT NULL,
  value                       TEXT NOT NULL,

  PRIMARY KEY (gene_expression_id),
  UNIQUE KEY gene_expression_idx(gene_id, tissue_id)



) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


CREATE TABLE tissue (
  tissue_id                   INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  name                        VARCHAR(255),
  description                 TEXT,

  PRIMARY KEY (tissue_id),
  UNIQUE KEY name_idx (name)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;
