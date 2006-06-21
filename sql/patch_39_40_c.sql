# patch_39_40_c
#
# title: xref version 
#
# description:
# Change version column in xref to be default 0, not allow null

ALTER TABLE xref CHANGE COLUMN version version VARCHAR(10) DEFAULT '0' NOT NULL;

UPDATE xref SET version = '0' WHERE version = NULL OR version = '' OR version = "NULL";

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_39_40_c.sql|xref_version');

