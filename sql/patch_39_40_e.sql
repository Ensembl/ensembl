# patch_39_40_e
#
# title: is_current not null
#
# description:
# this patch changes the is_current column in gene, transcript and exon to be
# 'NOT NULL', to prevent loading data where this property is not set

# change column 'is_current'

ALTER TABLE gene CHANGE is_current is_current BOOLEAN NOT NULL DEFAULT 1;
ALTER TABLE transcript CHANGE is_current is_current BOOLEAN NOT NULL DEFAULT 1;
ALTER TABLE exon CHANGE is_current is_current BOOLEAN NOT NULL DEFAULT 1;

# set is_current to 1
UPDATE gene set is_current = 1;
UPDATE transcript set is_current = 1;
UPDATE exon set is_current = 1;

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_39_40_e.sql|is_current_not_null');

