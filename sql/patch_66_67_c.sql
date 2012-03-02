# patch_66_67_c.sql
#
# Title: Adding intron_supporting_evidence table
#
# Description: Introns can be supported by an external feature. This gives a 
# weight to how much we believe the intron
# 

CREATE TABLE intron_supporting_evidence (
  intron_supporting_evidence_id INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  previous_exon_id INT(10) UNSIGNED NOT NULL,
  next_exon_id INT(10) UNSIGNED NOT NULL,
  hit_name VARCHAR(100) NOT NULL,
  score DECIMAL(10,3),
  score_type ENUM('NONE', 'DEPTH') DEFAULT 'NONE',
  
  PRIMARY KEY (intron_supporting_evidence_id),
  
  UNIQUE KEY (previous_exon_id, next_exon_id)
) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

# Patch identifier:
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_66_67_c.sql|adding_intron_supporting_evidence');
