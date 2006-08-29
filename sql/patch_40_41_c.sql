# patch_40_41_c
#
# title: xref priority
#
# description:
# Add a priority column to the xref table

ALTER TABLE xref ADD COLUMN priority INT DEFAULT 1 NOT NULL;

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_40_41_c.sql|xref_priority');