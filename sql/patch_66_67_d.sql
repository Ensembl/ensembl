# patch_66_67_d.sql
#
# Title: Adding new status type ANNOTATED
#
# Description: Including a new status type for genes and transcript
# 

ALTER TABLE gene MODIFY COLUMN status 
  ENUM('KNOWN', 'NOVEL', 'PUTATIVE', 'PREDICTED', 'KNOWN_BY_PROJECTION', 'UNKNOWN', 'ANNOTATED');
  
ALTER TABLE transcript MODIFY COLUMN status 
  ENUM('KNOWN', 'NOVEL', 'PUTATIVE', 'PREDICTED', 'KNOWN_BY_PROJECTION', 'UNKNOWN', 'ANNOTATED');

# Patch identifier:
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_66_67_d.sql|adding_gene_transcript_annotated');
