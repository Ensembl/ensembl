# patch_59_60_c.sql
#
# Title:
#   A patch to fix a couple of inconsistencies in the schema.
#
# Description:
#   A couple of fixes to the schema to do with inconsistencies found
#   during QC.

# Make the the 'seq_region_start' and 'seq_region_end' fields of the
# 'karyotype' table UNSIGNED (like they are everywhere else).
ALTER TABLE karyotype
  MODIFY COLUMN seq_region_start INT(10) UNSIGNED NOT NULL,
  MODIFY COLUMN seq_region_end   INT(10) UNSIGNED NOT NULL;

# Make 'seq_region.length' UNSIGNED (we do not like negative lengths).
ALTER TABLE seq_region
  MODIFY COLUMN length INT(10) UNSIGNED NOT NULL;

# Insert patch identifier.
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_59_60_c.sql|fix_inconsistencies');
