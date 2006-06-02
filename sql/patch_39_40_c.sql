# patch_39_40_c
#
# title: xref version 
#
# description:
# Change version column in xref to be default 0, not allow null

ALTER TABLE xref CHANGE COLUMN version version VARCHAR(10) DEFAULT '0' NOT NULL;

UPDATE xref SET version = '0' WHERE version = NULL OR version = '';

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_39_40_c.sql|xref_version');

























































# Similarly make all seq_region_start & seq_region_end columns INT(10) UNSIGNED
# Although these are not used as keys, they are used in many joins and having
# the same column type should make joins faster as MySQL will not have to do
# any type casting.








