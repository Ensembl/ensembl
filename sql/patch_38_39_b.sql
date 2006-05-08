# patch_38_39_b
#
# title: unique assembly
#
# description:
# this patch adds a unique index to the assembly table to prevent duplicate data

# create a new table with the unique index

CREATE TABLE assembly_new (

  asm_seq_region_id           INT UNSIGNED NOT NULL,
  cmp_seq_region_id           INT(10) UNSIGNED NOT NULL, 
  asm_start                   INT(10) NOT NULL,
  asm_end                     INT(10) NOT NULL,
  cmp_start                   INT(10) NOT NULL,
  cmp_end                     INT(10) NOT NULL,
  ori                         TINYINT  NOT NULL, 
  
  KEY (cmp_seq_region_id),
  KEY (asm_seq_region_id, asm_start),
  UNIQUE KEY all_idx (asm_seq_region_id, cmp_seq_region_id, asm_start, asm_end, cmp_start, cmp_end, ori)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;

# insert unique values into the new assembly table
INSERT IGNORE INTO assembly_new SELECT * FROM assembly;

# drop old assembly table and rename new one
DROP TABLE assembly;
ALTER TABLE assembly_new RENAME assembly;

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_38_39_b.sql|unique_assembly');

