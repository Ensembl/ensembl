# patch_41_42_a
#
# title: remove xref.priority
#
# description:
# Remove the priority column in xref

ALTER TABLE xref DROP COLUMN priority;

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_41_42_a.sql|remove_xref_priority');