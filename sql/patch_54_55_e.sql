# patch_54_55_e.sql
#
# title: Add column 'is_constitutive' to the exon table.
#
# description:
# The 'is_constitutive' column in the exon table will be set to true
# for exons that are constitutive.  This is done by a script in
# 'misc-scripts'.

ALTER TABLE exon ADD COLUMN is_constitutive BOOLEAN NOT NULL DEFAULT 0;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
VALUES (NULL, 'patch', 'patch_54_55_e.sql|add_is_constitutive_column');
